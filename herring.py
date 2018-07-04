import subprocess, os, copy, logging, itertools, sys, pickle
from Bio import SeqIO
logging.basicConfig(filename='GIfinder.log', level=logging.DEBUG, format='%(message)s')
sys.path.insert(0, '/home/asier/PycharmProjects/flex2')
import blastParser

#Windows stretcher path
#stretcherPath = 'C:\mEMBOSS\stretcher.exe'
stretcherPath = '/usr/bin/stretcher'
blastPath = '/opt/ncbi-blast-2.6.0+/bin/'
fastaList = ['WCFS1.fasta', 'LP2.fasta']
fastaListLengths = {}
refFasta = 'lactobacillus.core.fasta'

fileList = []
masterFamilyList = []
masterSeqList = []


class GapFamily():
    def __init__(self, gapList):
        self.parents = (gapList[0][0], gapList[0][4])
        self.gapList = []
        for gapItem in gapList:
            newGap = Gap(gapItem, self)
            self.gapList.append(newGap)

    def writeGapInfo(self, gapList, filename):
        with open(filename, 'w') as output:
            for i, gap in enumerate(gapList):
                gapStr = gap.buildInfoStr()
                output.write(gapStr)
                if i < len(gapList) - 1:
                    output.write('-----------------------------------\n\n')


    def createStrList(self):
        strList = []
        for i, gap in enumerate(self.gapList):
            strList.append(gap.name)
        return strList

    def findAdjacentGaps(self, gap):
        #sort gaps by seq1pos
        self.gapList.sort(key=lambda gap: gap.parentPos[0])
        #get the index for the target gap
        targetIndex = self.gapList.index(gap)
        prevGap = None
        nextGap = None
        logging.debug('target index: {}, gapListLen: {}'.format(targetIndex, len(self.gapList)))
        if targetIndex != len(self.gapList)-1:
            nextGap = targetIndex + 1
        if targetIndex != 0:
            prevGap = targetIndex - 1

        if prevGap is not None:
            prevGap = self.gapList[prevGap]
        if nextGap is not None:
            nextGap = self.gapList[nextGap]

        return (prevGap, nextGap)

    def equalize(self):
        refName = getRefFastaName(refFasta)
        if self.parents[0] == refName:
            parents1 = self.parents[0]
            self.parents = (refName, parents1)

        for gap in self.gapList:
            if gap.parentSeqs[0] != self.parents[0]:
                newParent = gap.corePos
                newCore = gap.parentPos

                gap.parentSeqs = self.parents
                gap.corePos = newCore
                gap.parentPos = newParent

    def __iter__(self):
        return iter(self.gapList)

    def closeGaps(self, rejectedGaps = ()):
        #recursive algorithm,looks completely insane but it might work?
        self.gapList.sort(key=lambda gap: gap.parent1pos[0])
        for gap in self.gapList:
            if gap.isGapOrphan() == [True, 'left'] and gap not in rejectedGaps:
                targetIndex = self.gapList.index(gap)
                for i in range(targetIndex + 1, len(self.gapList)):
                    if self.gapList[i].isGapOrphan[0] == True:
                        for borderGap in gap.leftMatch:
                            if borderGap in self.gapList[i].rightMatch:
                                newGapList = [gap.parents[0], gap.parent1pos[0], self.gapList[i].parent1pos[1], None,
                                              gap.parents[1], gap.parent2pos[0], self.gapList[i].parent2pos[1], None,
                                              gap.type]
                                newGap = Gap(newGapList, self)
                                self.gapList.append(newGap)
                                self.gapList.remove(gap)
                                self.gapList.remove(self.gapList[i])
                                self.closeGaps(rejectedGaps)
                rejectedGaps.append(gap)
                self.closeGaps(rejectedGaps)

    def closeGaps2(self):
        self.gapList.sort(key=lambda gap: gap.parent1pos[0])
        for gap in self.gapList:
            if gap.isGapOrphanLeft() > 0:
                borderList = gap.isGapOrphanLeft()
                targetIndex = self.gapList.index(gap)
                for i in range(targetIndex + 1, len(self.gapList)):
                    pass

    def findOrphans(self, type = ('left', 'right', 'both')):
        orphanList = []
        for gap in self.gapList:
            orphanInfo = gap.isGapOrphan()
            if orphanInfo[0] == True and orphanInfo[1] in type:
                orphanList.append(gap)
        return orphanList

class Gap():
    def __init__(self, gapList, family):
        self.family = family
        self.parentSeqs = (gapList[0], gapList[4])
        self.corePos = (gapList[1], gapList[2])
        self.parentPos = (gapList[5], gapList[6])
        self.type = gapList[8]
        self.leftBorder = None
        self.rightBorder = None
        self.leftMatch = []
        self.rightMatch = []
        no = str(len(self.family.gapList))
        self.name = str(str(self.parentSeqs[0]) + '-' + str(self.parentSeqs[1]) + '_' + no)

    def buildInfoStr(self):
        diagStr = ''
        diagStr += ('Gap Name: {}\n'.format(self.name))
        diagStr += ('PARENT INFO:\n')
        diagStr += ('Parent 1: {}\n\t{} - {}\n'.format(self.parentSeqs[0], self.corePos[0], self.corePos[1]))
        diagStr += ('Parent 2: {}\n\t{} - {}\n'.format(self.parentSeqs[1], self.parentPos[0], self.parentPos[1]))
        diagStr += ('EDGES\nLEFT BORDER:\n')
        for match in self.leftMatch:
            diagStr += ('\t{}\n'.format(match.name))
        diagStr += ('\nRIGHT BORDER:\n')
        for match in self.rightMatch:
            diagStr += ('\t{}\n'.format(match.name))
        return diagStr

    def getParentSeq(self, parent=1):
        gapName = ''
        gapSeq = None
        # get sequence from the fasta file
        for filePair in fileList:
            if self.parentSeqs[parent] in filePair:
                gapName = filePair[parent]
        for record in SeqIO.parse(gapName, 'fasta'):
            if record.name == self.parentSeqs[1]:
                gapSeq = record.seq
        return gapSeq

    def getBorderSequence(self, size):
        gapName = ''
        gapSeq = self.getParentSeq()
        gapSeqLeft = (0, None)
        gapSeqRight = (0, None)

        #get sequence from the fasta file
        '''
        for filePair in fileList:
            if self.parents[0] in filePair:
                gapName = filePair[1]
        for record in SeqIO.parse(gapName, 'fasta'):
            if record.name == self.parents[0]:
                gapSeq = record.seq
        '''

        #Check if we can pick a sequence of the required size without taking a gap by mistake:
        adjGaps = self.family.findAdjacentGaps(self)

        #Let's start with the left border
        if adjGaps[0] is not None:
            leftGapPos = adjGaps[0].parentPos[1]
            if (self.parentPos[0] - leftGapPos) < size:
                logging.debug('there is not enough space for a full left border! size: {}, space between gaps: {}'.format(size, self.parentPos[0] - leftGapPos))
                gapSeqLeft = ((self.parentPos[0] - leftGapPos), gapSeq[leftGapPos:self.parentPos[0]])
            else:
                gapSeqLeft = (size, gapSeq[(self.parentPos[0] - size):self.parentPos[0]])


        #right border is next
        if adjGaps[1] is not None:
            rightGapPos = adjGaps[1].parentPos[0]
            if (rightGapPos - self.parentPos[1]) < size:
                logging.debug('there is not enough space for a full right border! size: {}, space between gaps: {}'.format(size,rightGapPos - self.parentPos[1]))
                gapSeqRight = ((rightGapPos - self.parentPos[1]), gapSeq[self.parentPos[1]:rightGapPos])
            else:
                gapSeqRight = (size, gapSeq[self.parentPos[1]:self.parentPos[1] + size])

        #then return both borders
        self.leftBorder = gapSeqLeft
        self.rightBorder = gapSeqRight
        return [gapSeqLeft, gapSeqRight]

    def addBorderGapRef(self, gap, type):
        if type == 'left':
            self.leftMatch.append(gap)
        elif type == 'right':
            self.rightMatch.append(gap)

    def isGapOrphan(self):
        if len(self.leftMatch) + len(self.rightMatch) == 0:
            return [True, 'both']
        elif len(self.leftMatch) == 0 and len(self.rightMatch) > 0:
            return [True, 'right']
        elif len(self.leftMatch) > 0 and len(self.rightMatch) == 0:
            return [True, 'left']
        else:
            borderToCheck = [self.leftMatch, self.rightMatch]
            if len(self.leftMatch) < len(self.rightMatch):
                borderToCheck = [self.rightMatch, self.leftMatch]
            for borderGap in borderToCheck[0]:
                if borderGap in borderToCheck[1]:
                    return [False, None]
            return [True, 'both']

def saveMasterFamilyList(masterFamilyList, filename):
    with open(filename, 'wb') as output:
        pickle.dump(masterFamilyList, output, pickle.HIGHEST_PROTOCOL)

def loadMasterFamilyList(filename):
    with open(filename, 'rb') as input:
        return pickle.load(input)

def alignSeqsByStretcher(seq1, seq2, gapOpen = 1000, gapExtend = 1000):
    logging.info('Aligning 2 Sequences')
    outfile = 'stretcher.aln.temp'
    logging.info('Running stretcher...')
    subprocess.call([stretcherPath, '-asequence', seq1, '-bsequence', seq2, '-outfile', outfile, '-gapopen', str(gapOpen), '-gapextend',
          str(gapExtend)])
    logging.info('Stretcher finished, extracting alignment stats from alignment file')
    statsDict2 = {}
    with open(outfile) as outHandle:
        # 23 is identity, 24 is similarity, 25 is gaps, 26 is score
        statsDict = {'23': 'Identity', '24': 'Similarity', '25': 'Gaps'}

        for i, line in enumerate(outHandle):
            if i < 30:
                if str(i) in statsDict.keys():
                    lineSplit = line.split('(')[-1]
                    lineSplit = float(lineSplit[0:-3])
                    statsDict2[statsDict[str(i)]] = lineSplit
        os.remove(outfile)
    logging.info('Alignment Stats are as follows: {}'.format(statsDict2))
    return statsDict2

def findGaps(blastFamily, thresholdMin, thresholdMax, masterFamilyList):
    logging.info('Looking for gaps...')
    logging.info('ThresholdMin is {}, ThresholdMax is {}'.format(thresholdMin, thresholdMax))

    blastFamily._equalize()
    #family.rearrangeBlastList()
    blastFamily.sortHits()

    logging.info('N of total Blasts: {}'.format(len(blastFamily.blastList)))
    blastList = blastFamily.blastList
    gapList = []

    for i in range(0, len(blastList) - 1):
        logging.debug('Iteration n {}'.format(i))
        total = 0
        #Get the 2 blasts to compare
        fstBlast = copy.copy(blastList[i])
        scdBlast = copy.copy(blastList[i + 1])

        #Check if the blasts are reversed, and if they are fix them
        if fstBlast.seq1pos[1] < fstBlast.seq1pos[0] or fstBlast.seq1pos[1] < fstBlast.seq1pos[0]:
            fstBlast.seq1pos = (fstBlast.seq1pos[1], fstBlast.seq1pos[0])
        if scdBlast.seq1pos[1] < scdBlast.seq1pos[0] or scdBlast.seq1pos[1] < scdBlast.seq1pos[0]:
            scdBlast.seq1pos = (scdBlast.seq1pos[1], scdBlast.seq1pos[0])

        #calculate distance between both blasts in both sequences
        logging.info('blast1pos: {} - {}'.format(fstBlast.seq1pos[0], fstBlast.seq1pos[1]))
        logging.info('blast2pos: {} - {}'.format(scdBlast.seq1pos[0], scdBlast.seq1pos[1]))
        seq1dtce = (scdBlast.seq1pos[0] - fstBlast.seq1pos[1])
        seq2dtce = (scdBlast.seq2pos[0] - fstBlast.seq2pos[1])
        logging.debug('Distances between blasts: seq1 {}, seq2 {}'.format(seq1dtce, seq2dtce))

        #Check if the distances in both sequences meet the required threshold
        if seq1dtce > thresholdMin and seq1dtce < thresholdMax:
            total += 1
        if seq2dtce > thresholdMin and seq2dtce < thresholdMax:
            total += 1

        #then assign a category depending on the distances
        #gapList is: parent1, seq1pos1, seq1pos2, seq1dtce, parent2, seq2pos1, seq2pos2, seq2dtce, type
        if total == 0:
            logging.debug('Total is {}, which means there is no gap'.format(total))
            pass
        elif total == 1:
            logging.debug('Total is {}, which means there is a gap, probably an indel'.format(total))
            gapList.append([fstBlast.parents[0], fstBlast.seq1pos[1], scdBlast.seq1pos[0],seq1dtce, fstBlast.parents[1], fstBlast.seq2pos[1], scdBlast.seq2pos[0],seq2dtce, 'indel'])
        elif total == 2:
            logging.debug('Total is {}, which means there is a gap, probably a gap'.format(total))
            gapList.append([fstBlast.parents[0], fstBlast.seq1pos[1], scdBlast.seq1pos[0],seq1dtce, fstBlast.parents[1], fstBlast.seq2pos[1], scdBlast.seq2pos[0],seq2dtce, 'gap'])
        logging.debug('Sotring gap as gapList {}'.format([fstBlast.parents[0], fstBlast.seq1pos[1], scdBlast.seq1pos[0],seq1dtce, fstBlast.parents[1], fstBlast.seq2pos[1], scdBlast.seq2pos[0],seq2dtce, 'type']))

    #Sort the list by seq1pos1
    gapList.sort(key=lambda x : x[1])
    #create gap family, then add it to the master list
    print(gapList)
    if len(gapList) > 0:
        newGapFamily = GapFamily(gapList)
        masterFamilyList.append(newGapFamily)

def getRefFastaName(refFasta):
    return SeqIO.read(refFasta, 'fasta').name

def createCombinedFasta(fastaList):
    logging.info('Combining fastas')
    with open('combinedFasta.temp', 'w') as output:
        for fasta in fastaList:
            for record in SeqIO.parse(fasta, 'fasta'):
                fileList.append((str(record.name), fasta))
                if record.name not in fastaListLengths.keys():
                    fastaListLengths[record.name] = len(record.seq)
                output.write('>' + str(record.name) + '\n')
                output.write(str(record.seq) + '\n')
    logging.info('Fastas combined')

def removeBlastDBFiles():
    logging.info('cleaning Blast DB files...')
    dirFiles = os.listdir(os.curdir)
    removedFileList = []
    removedFileCounter = 0
    for file in dirFiles:
        if file.split('.')[0] == 'dbTemp':
            removedFileCounter += 1
            removedFileList.append(file)
            os.remove(file)
    logging.info('{} files removed: {}'.format(removedFileCounter, removedFileList))

def performBlastAgainstCore(fastaList, refFasta, minThres, maxThres, masterFamilyList):
    for fasta in fastaList:
        createCombinedFasta([fasta, refFasta])
        subprocess.call([blastPath + 'makeblastdb', '-in', 'combinedFasta.temp', '-out', 'dbTemp', '-dbtype', 'nucl'])
        subprocess.call(
            [blastPath + 'blastn', '-query', 'combinedFasta.temp', '-db', 'dbTemp', '-out', 'blastSeqs.blastn',
             '-num_threads', '4', '-outfmt', '6'])
        os.remove('combinedFasta.temp')
        removeBlastDBFiles()

        blastHits = blastParser.parseBlastFile('blastSeqs.blastn', minIdentity=90, minAln=1500)
        blastFamilies = blastParser.groupHits(blastHits)
        for family in blastFamilies:
            family._equalize()
            family.addParentLengths(fastaListLengths)
            family.removeOwnHits()
            family.removeInternalHits()
            family.removeSamePosHits()
            family.removeStrangeHits()
            findGaps(family, minThres, maxThres, masterFamilyList)

def equalizeBorderSizes(gap1borders, gap2borders):

    # remember: getBorderSequence returns a tuple: ((leftBorderSize, leftBorderSeq), (rightBorderSize, rightBorderSeq))
    # borders might have different sizes, so we have to equalize the border size
    # start with the left border:
    if (gap1borders[0][0] + gap2borders[0][0]) > 0:
        if gap1borders[0][0] != gap2borders[0][0]:
            logging.info('left side borders are not of the same size! gap1: {}, gap2: {}'.format(gap1borders[0][0], gap2borders[0][0]))

            if gap1borders[0][0] > gap2borders[0][0]:
                #we calculate the difference with the total bordar size because it has the same result
                diff = gap1borders[0][0] - gap2borders[0][0]
                gap1borders[0] = (gap1borders[0][0] - diff, gap1borders[0][1][diff:])
                logging.info('Fixed gap1 border: now gap1size is {} and gap2size is {}'.format(gap1borders[0][0], gap2borders[0][0]))
            elif gap2borders[0][0] > gap1borders[0][0]:
                diff = gap2borders[0][0] - gap1borders[0][0]
                gap2borders[0] = (gap2borders[0][0] - diff, gap2borders[0][1][diff:])
                logging.info(
                    'Fixed gap2 border: now gap1size is {} and gap2size is {}'.format(gap1borders[0][0], gap2borders[0][0]))
            logging.debug('gap1leftBorder: {},  gap2leftBorder: {}'.format(len(gap1borders[0][1]), len(gap2borders[0][1])))


    #right border:
    elif (gap1borders[1][0] + gap2borders[1][0]) > 0:
        if gap1borders[1][0] != gap2borders[1][0]:
            logging.debug('right side borders are not of the same size! gap1: {}, gap2: {}'.format(gap1borders[1][0],
                                                                                                    gap2borders[1][0]))
            if gap1borders[1][0] > gap2borders[1][0]:
                diff = gap1borders[1][0] - gap2borders[1][0]
                gap1borders[1] = (gap1borders[1][0] - diff, gap1borders[1][1][:len(gap1borders[1][1]) - diff])
                logging.info(
                    'Fixed gap1 border: now gap1size is {} and gap2size is {}'.format(gap1borders[1][0], gap2borders[1][0]))
            elif gap2borders[1][0] > gap1borders[1][0]:
                diff = gap2borders[1][0] - gap1borders[1][0]
                gap2borders[1] = (gap2borders[1][0] - diff, gap2borders[1][1][:len(gap2borders[1][1]) - diff])
                logging.info(
                    'Fixed gap2 border: now gap1size is {} and gap2size is {}'.format(gap1borders[1][0], gap2borders[1][0]))
            logging.debug('gap1rightBorder: {}, gap2rightBorder: {}'.format(len(gap2borders[1][1]), len(gap2borders[1][1])))

        else:
            logging.info('No need to fix borders: gap1size is {} and gap2size is {}'.format(gap1borders[1][0], gap2borders[1][0]))

def compareTwoGaps(gap1, gap2, size, storeResults = True):
    # Get the border sequences for both gaps
    gap1borders = gap1.getBorderSequence(size)
    gap2borders = gap2.getBorderSequence(size)
    if gap1borders[0][1] != None and gap2borders[0][1] != None:
        equalizeBorderSizes(gap1borders, gap2borders)



    #call stretcher and get the alignment. first, prepare the fasta files
    tempFastaFiles = ['gap1left.temp.fasta', 'gap1right.temp.fasta', 'gap2left.temp.fasta', 'gap2right.temp.fasta']
    rightalnResults = None
    leftalnResults = None
    logging.info('comparing left border')
    if gap1borders[0][1] != None and gap2borders[0][1] != None:
        with open(tempFastaFiles[0], 'w') as gap1left:
            gap1left.write('>gap1left\n')
            gap1left.write(str(gap1borders[0][1]))
        with open(tempFastaFiles[2], 'w') as gap2left:
            gap2left.write('>gap2left\n')
            gap2left.write(str(gap2borders[0][1]))
        leftalnResults = alignSeqsByStretcher('gap1left.temp.fasta', 'gap2left.temp.fasta')
    else:
        logging.info('No sequences for the left border')
        leftalnResults = {'Identity':-1, 'Gaps':0}

    logging.info('comparing right border')
    if gap1borders[1][1] != None and gap2borders[1][1] != None:
        with open(tempFastaFiles[1], 'w') as gap1right:
            gap1right.write('>gap1right\n')
            gap1right.write(str(gap1borders[1][1]))
        with open(tempFastaFiles[3], 'w') as gap2right:
            gap2right.write('>gap2right\n')
            gap2right.write(str(gap2borders[1][1]))
        rightalnResults = alignSeqsByStretcher('gap1right.temp.fasta', 'gap2right.temp.fasta')
    else:
        logging.info('No sequences for the right border')
        rightalnResults = {'Identity': -1, 'Gaps': 0}


    #check if they are similar. If they are, then add borders
    response = [False, False]
    #todo - check if reverse complement sequences fuck with identity(use blast then)
    if leftalnResults['Identity'] > 95 and leftalnResults['Gaps'] < 1:
        logging.debug('left sides match')
        response[0] = True
        if storeResults == True:
            gap1.addBorderGapRef(gap2, 'left')
            gap2.addBorderGapRef(gap1, 'left')

    if rightalnResults['Identity'] > 95 and rightalnResults['Gaps'] < 1:
        logging.debug('right sides match')
        response[0] = True
        if storeResults == True:
            gap1.addBorderGapRef(gap2, 'right')
            gap2.addBorderGapRef(gap1, 'right')

    #remove all temp fasta files
    for file in tempFastaFiles:
        if file in os.listdir(os.curdir):
            os.remove(file)

    return response

def compareGapFamilies(gapFamily1, gapFamily2, size):
    # we want all combinations, not permutations (permutations also include the order of the elements - so AB & BA)
    # combList = itertools.combinations(gapFamily1.gapList + gapFamily2.gapList, 2)

    # Actually we don't want all combinations either, we just want all combinations between 2 sets (e.g. if list1 = 25 and list2 = 15,
    # there would be 780 combinations but only 300 combinations between both sets -> Try list comprehensions instead, and build an iterator
    # if it becomes too much of a memory burden (it will)
    totalNoOfCombinations = len(gapFamily1.gapList) * len(gapFamily2.gapList)
    currComb = 0
    combList = iter([(x, y) for x in gapFamily1.gapList for y in gapFamily2.gapList])

    for gapPair in combList:
        currComb += 1
        print('{} / {}'.format(currComb, totalNoOfCombinations))
        if gapPair[0].family != gapPair[1].family:
            logging.debug('comparing gap pair {} - {}'.format(gapPair[0].name, gapPair[1].name))
            compareTwoGaps(gapPair[0], gapPair[1], size)


if __name__ == '__main__':
    masterFamilyList = []
    '''
    blastFamilies = performBlastAgainstCore(fastaList, refFasta, 2000, 1000000, masterFamilyList)
    for family in masterFamilyList:
        family.equalize()

    #Perform comparisons
    familyCombList = itertools.combinations(masterFamilyList, 2)
    for familyPair in familyCombList:
        compareGapFamilies(familyPair[0], familyPair[1], 5000)
    '''
    masterFamilyList = loadMasterFamilyList('masterListTest.pk1')
    for family in masterFamilyList:
        print(len(family.gapList))
        family.writeGapInfo(family.findOrphans(type=['left']), 'diagnostic-Left' + '-' + str(family.parents[0]) + '.txt')




