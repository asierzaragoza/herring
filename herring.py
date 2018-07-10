#todo - right now the fasta has to have the same name as its header to work. Not a critical error, but a pain in the ass

import subprocess, os, copy, logging, itertools, sys, pickle, time
from Bio import SeqIO

#Log setup so we can have more than 1 log at the same time (shamelessly stolen from stackOverflow)
def setup_logger(logger_name, log_file, level=logging.DEBUG):
    l = logging.getLogger(logger_name)
    formatter = logging.Formatter('%(message)s')
    fileHandler = logging.FileHandler(log_file, mode='a')
    fileHandler.setFormatter(formatter)
    #streamHandler = logging.StreamHandler()
    #streamHandler.setFormatter(formatter)

    l.setLevel(level)
    l.addHandler(fileHandler)
    #l.addHandler(streamHandler)

setup_logger('herringMain', 'herring.main.log')
loggerMain = logging.getLogger('herringMain')


sys.path.insert(0, '/home/asier/PycharmProjects/flex2')
import blastParser

#Windows stretcher path
#stretcherPath = 'C:\mEMBOSS\stretcher.exe'
stretcherPath = '/usr/bin/stretcher'
#blastPath = '/opt/ncbi-blast-2.6.0+/bin/'
blastPath = '/home/rohit/asier/programs/ncbi-blast-2.7.1+/bin/'
prodigalPath = '/home/rohit/asier/programs/prodigal/'
fastaList = ['WCFS1.fasta', 'LP2.fasta', 'B21.fasta']
fastaListLengths = {}
refFasta = 'lactobacillus.coreWCFS1.fasta'

fileList = []
masterFamilyList = []



class GapFamily():
    def __init__(self, gapList):
        self.parents = (gapList[0][0], gapList[0][4])
        self.gapList = []
        for gapItem in gapList:
            newGap = Gap(gapItem, self)
            self.gapList.append(newGap)

    def writeGapInfo(self, gapList, filename):
        with open(filename, 'w') as output:
            output.write('NUMBER OF GAPS: {}'.format(len(gapList)))
            for i, gap in enumerate(gapList):

                gapStr = gap.buildInfoStr()
                output.write(gapStr)
                if i < len(gapList) - 1:
                    output.write('-----------------------------------\n\n')

    def writeGapInfo2(self, gapList):
        pass

    def getParentSeq(self, parent=0):
        gapName = ''
        gapSeq = None
        # get sequence from the fasta file
        for filePair in fileList:
            if self.parents[parent] in filePair:
                gapName = filePair[parent]
        for record in SeqIO.parse(gapName, 'fasta'):
            if record.name == self.parents[1]:
                gapSeq = record.seq
        return gapSeq

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
        loggerMain.debug('target index: {}, gapListLen: {}'.format(targetIndex, len(self.gapList)))
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
        loggerMain.info('parents[0] = {}, refName = {}'.format(self.parents[0], refName))
        if self.parents[0] == refName:
            parents1 = self.parents[0]
            self.parents = (refName, parents1)
            loggerMain.info('fixed: parent[0] (parent) is {} and parent[1] (core) is {}'.format(self.parents[0], self.parents[0]))
        for gap in self.gapList:
            loggerMain.info('gap.parents[0] = {}, family.parents[0] is {}'.format(gap.parents[0], self.parents[0]))
            if gap.parents[0] != self.parents[0]:

                newParent = copy.copy(gap.corePos)
                newCore = copy.copy(gap.parentPos)

                gap.parents = self.parents
                gap.corePos = newCore
                gap.parentPos = newParent
            loggerMain.info('parents = {} , corePos = {} , parentPos = {}'.format(gap.parents, gap.corePos, gap.parentPos))

    def __iter__(self):
        return iter(self.gapList)

    def getGapsByType(self, type='analyze'):
        targetList = []
        for gap in self.gapList:
            if gap.type == type:
                targetList.append(gap)
        return targetList

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
        self.parents = (gapList[0], gapList[4])
        self.corePos = (gapList[5], gapList[6])
        self.parentPos = (gapList[1], gapList[2])
        self.type = gapList[8]
        self.leftBorder = None
        self.rightBorder = None
        self.leftMatch = []
        self.rightMatch = []
        no = str(len(self.family.gapList))
        self.name = str(str(self.parents[0]) + '-' + str(self.parents[1]) + '_' + no)

    def buildInfoStr(self):
        diagStr = ''
        diagStr += ('Gap Name: {}\n'.format(self.name))
        diagStr += ('PARENT INFO:\n')
        diagStr += ('Parent 1: {}\n\t{} - {}\n'.format(self.parents[1], self.corePos[0], self.corePos[1]))
        diagStr += ('Parent 2: {}\n\t{} - {}\n'.format(self.parents[0], self.parentPos[0], self.parentPos[1]))
        diagStr += ('EDGES\nLEFT BORDER: {}\n'.format(len(self.leftMatch)))
        for match in self.leftMatch:
            diagStr += ('\t{}\n'.format(match.name))
        diagStr += ('\nRIGHT BORDER: {}\n'.format(len(self.rightMatch)))
        for match in self.rightMatch:
            diagStr += ('\t{}\n'.format(match.name))
        return diagStr

    def getParentSeq(self, parent=0):
        gapName = ''
        gapSeq = None
        # get sequence from the fasta file
        for filePair in fileList:
            if self.parents[parent] in filePair:
                gapName = filePair[parent]
        for record in SeqIO.parse(gapName + '.fasta', 'fasta'):
            if record.name == gapName:
                gapSeq = record.seq
        return gapSeq

    def getBorderSequence(self, size):
        loggerMain.info('finding border sequences for {} with size {}'.format(self.name, size))
        loggerMain.info('parent sequence used is {}'.format(self.parents[0]))
        gapSeq = self.getParentSeq()
        gapSeqLeft = (0, None)
        gapSeqRight = (0, None)



        #Check if we can pick a sequence of the required size without taking a gap by mistake:
        adjGaps = self.family.findAdjacentGaps(self)

        #Let's start with the left border
        if adjGaps[0] is not None:
            leftGapPos = adjGaps[0].parentPos[1]
            if (self.parentPos[0] - leftGapPos) < size:
                loggerMain.debug('there is not enough space for a full left border! size: {}, space between gaps: {}'.format(size, self.parentPos[0] - leftGapPos))
                gapSeqLeft = ((self.parentPos[0] - leftGapPos), gapSeq[leftGapPos:self.parentPos[0]])
            else:
                gapSeqLeft = (size, gapSeq[(self.parentPos[0] - size):self.parentPos[0]])
        else:
            if (self.parentPos[0] - size) > 0:
                gapSeqLeft = (size, gapSeq[(self.parentPos[0] - size):self.parentPos[0]])
            else:
                gapSeqLeft = (size, gapSeq[:self.parentPos[0]])


        #right border is next
        if adjGaps[1] is not None:
            rightGapPos = adjGaps[1].parentPos[0]
            if (rightGapPos - self.parentPos[1]) < size:
                loggerMain.debug('there is not enough space for a full right border! size: {}, space between gaps: {}'.format(size,rightGapPos - self.parentPos[1]))
                gapSeqRight = ((rightGapPos - self.parentPos[1]), gapSeq[self.parentPos[1]:rightGapPos])
            else:
                gapSeqRight = (size, gapSeq[self.parentPos[1]:self.parentPos[1] + size])
        else:
            if (self.parentPos[1] + size) < len(gapSeq):
                gapSeqRight = (size, gapSeq[self.parentPos[1]:self.parentPos[1] + size])
            else:
                gapSeqRight = (size, gapSeq[self.parentPos[1]:])

        #then return both borders
        self.leftBorder = gapSeqLeft
        self.rightBorder = gapSeqRight
        loggerMain.info('borders for {} found: left border is {}, while right border is {}'.format(self.name,[(self.parentPos[0] - self.leftBorder[0]), self.parentPos[0]], [self.parentPos[1], self.parentPos[1] + self.rightBorder[0]] ))
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


def saveHerringData(filename):
    with open(filename, 'wb') as output:
        pickle.dump(masterFamilyList, output, pickle.HIGHEST_PROTOCOL)
        pickle.dump(fileList, output, pickle.HIGHEST_PROTOCOL)

def loadHerringData(filename):
    with open(filename, 'rb') as input:
        familylist = pickle.load(input)
        filelist = pickle.load(input)
    return(familylist, filelist)

def alignSeqsByStretcher(seq1, seq2, gapOpen = 1000, gapExtend = 1000):
    loggerMain.info('Aligning 2 Sequences by Stretcher')
    outfile = './temp/stretcher.aln.temp'
    loggerMain.info('Running stretcher...')
    subprocess.call([stretcherPath, '-asequence', seq1, '-bsequence', seq2, '-outfile', outfile, '-gapopen', str(gapOpen), '-gapextend',
          str(gapExtend)])
    loggerMain.info('Stretcher finished, extracting alignment stats from alignment file')
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
    loggerMain.info('Alignment Stats are as follows: {}'.format(statsDict2))
    return statsDict2

def alignSeqsByBlast(seq1, seq2, borderSize):
    loggerMain.info('Aligning 2 Sequences by Blast')
    outfile = './temp/blastResults.temp.blastn'
    loggerMain.info('Running blastn...')
    subprocess.call([blastPath + 'blastn', '-query', seq1, '-subject', seq2, '-out', outfile, '-outfmt', '6'])
    statsDict = {'Matches%': 0, 'AlnLen': 0, 'Gaps': 0, 'Mismatches': 0}
    if os.stat(outfile).st_size > 0:
        loggerMain.info('There was a match!')
        with open(outfile, 'r') as blastHandle:
            statsDict = {'Matches%':0 , 'AlnLen': 0, 'Gaps': 0, 'Mismatches' : 0}
            #remember: order is queryID / subjectID / percentage of identical matches / alignment length / n of mismatches / n of gaps / start in query
            # end in query / start in subject / end in subject / E-value / bitscore
            blastLineInfo = blastHandle.readline().split('\t')
            statsDict['Matches%'] = float(blastLineInfo[2])
            statsDict['AlnLen'] = float(int(blastLineInfo[3]) / borderSize)
            statsDict['Mismatches'] = int(blastLineInfo[4])
            statsDict['Gaps'] = int(blastLineInfo[5])

    else:
        loggerMain.info('No match found')
    os.remove(outfile)
    return statsDict

def findAllGaps(blastFamily, thresholdMin, thresholdAnalysis, masterFamilyList):
    loggerMain.info('Looking for gaps...')
    loggerMain.info('ThresholdMin is {}, ThresholdAnalysis is {}'.format(thresholdMin, thresholdAnalysis))

    blastFamily._equalize()
    blastFamily.sortHits()

    loggerMain.info('N of total Blasts: {}'.format(len(blastFamily.blastList)))
    blastList = blastFamily.blastList
    gapList = []

    for i in range(0, len(blastList) - 1):
        loggerMain.debug('Iteration n {}'.format(i))

        # Get the 2 blasts to compare
        fstBlast = copy.copy(blastList[i])
        scdBlast = copy.copy(blastList[i + 1])

        # Check if the blasts are reversed, and if they are fix them
        if fstBlast.seq1pos[1] < fstBlast.seq1pos[0] or fstBlast.seq1pos[1] < fstBlast.seq1pos[0]:
            fstBlast.seq1pos = (fstBlast.seq1pos[1], fstBlast.seq1pos[0])
        if scdBlast.seq1pos[1] < scdBlast.seq1pos[0] or scdBlast.seq1pos[1] < scdBlast.seq1pos[0]:
            scdBlast.seq1pos = (scdBlast.seq1pos[1], scdBlast.seq1pos[0])

        # calculate distance between both blasts in both sequences
        loggerMain.info('blast1pos: {} - {}'.format(fstBlast.seq1pos[0], fstBlast.seq1pos[1]))
        loggerMain.info('blast2pos: {} - {}'.format(scdBlast.seq1pos[0], scdBlast.seq1pos[1]))
        seq1dtce = (scdBlast.seq1pos[0] - fstBlast.seq1pos[1])
        seq2dtce = (scdBlast.seq2pos[0] - fstBlast.seq2pos[1])
        loggerMain.debug('Distances between blasts: seq1 {}, seq2 {}'.format(seq1dtce, seq2dtce))
        total = 0
        # Check if the distances in both sequences meet the required threshold
        if seq1dtce > thresholdMin:
            total += 1
        if seq2dtce > thresholdMin:
            total += 1

        # then assign a category depending on the distances
        # gapList is: parent1, seq1pos1, seq1pos2, seq1dtce, parent2, seq2pos1, seq2pos2, seq2dtce, type
        if total == 0:
            loggerMain.debug('Total is {}, which means there is no gap'.format(total))
            pass
        elif total > 0:
            loggerMain.debug('Total is {}, which means there is a gap'.format(total))
            type = 'ignore'

            if seq1dtce + seq2dtce > thresholdAnalysis:
                type = 'analyze'

            gapList.append(
                [fstBlast.parents[0], fstBlast.seq1pos[1], scdBlast.seq1pos[0], seq1dtce, fstBlast.parents[1],
                 fstBlast.seq2pos[1], scdBlast.seq2pos[0], seq2dtce, type])
            loggerMain.debug('Sotring gap as gapList {}'.format(
                [fstBlast.parents[0], fstBlast.seq1pos[1], scdBlast.seq1pos[0], seq1dtce, fstBlast.parents[1],
                 fstBlast.seq2pos[1], scdBlast.seq2pos[0], seq2dtce, type]))

    # Sort the list by seq1pos1
    gapList.sort(key=lambda x: x[1])
    # create gap family, then add it to the master list
    if len(gapList) > 0:
        newGapFamily = GapFamily(gapList)
        masterFamilyList.append(newGapFamily)

def findGaps(blastFamily, thresholdMin, thresholdMax, masterFamilyList):
    loggerMain.info('Looking for gaps...')
    loggerMain.info('ThresholdMin is {}, ThresholdMax is {}'.format(thresholdMin, thresholdMax))

    blastFamily._equalize()
    #family.rearrangeBlastList()
    blastFamily.sortHits()

    loggerMain.info('N of total Blasts: {}'.format(len(blastFamily.blastList)))
    blastList = blastFamily.blastList
    gapList = []

    for i in range(0, len(blastList) - 1):
        loggerMain.debug('Iteration n {}'.format(i))
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
        loggerMain.info('blast1pos: {} - {}'.format(fstBlast.seq1pos[0], fstBlast.seq1pos[1]))
        loggerMain.info('blast2pos: {} - {}'.format(scdBlast.seq1pos[0], scdBlast.seq1pos[1]))
        seq1dtce = (scdBlast.seq1pos[0] - fstBlast.seq1pos[1])
        seq2dtce = (scdBlast.seq2pos[0] - fstBlast.seq2pos[1])
        loggerMain.debug('Distances between blasts: seq1 {}, seq2 {}'.format(seq1dtce, seq2dtce))

        #Check if the distances in both sequences meet the required threshold
        if seq1dtce > thresholdMin and seq1dtce < thresholdMax:
            total += 1
        if seq2dtce > thresholdMin and seq2dtce < thresholdMax:
            total += 1

        #then assign a category depending on the distances
        #gapList is: parent1, seq1pos1, seq1pos2, seq1dtce, parent2, seq2pos1, seq2pos2, seq2dtce, type
        if total == 0:
            loggerMain.debug('Total is {}, which means there is no gap'.format(total))
            pass
        elif total == 1:
            loggerMain.debug('Total is {}, which means there is a gap, probably an indel'.format(total))
            gapList.append([fstBlast.parents[0], fstBlast.seq1pos[1], scdBlast.seq1pos[0],seq1dtce, fstBlast.parents[1], fstBlast.seq2pos[1], scdBlast.seq2pos[0],seq2dtce, 'indel'])
        elif total == 2:
            loggerMain.debug('Total is {}, which means there is a gap, probably a gap'.format(total))
            gapList.append([fstBlast.parents[0], fstBlast.seq1pos[1], scdBlast.seq1pos[0],seq1dtce, fstBlast.parents[1], fstBlast.seq2pos[1], scdBlast.seq2pos[0],seq2dtce, 'gap'])
        loggerMain.debug('Sotring gap as gapList {}'.format([fstBlast.parents[0], fstBlast.seq1pos[1], scdBlast.seq1pos[0],seq1dtce, fstBlast.parents[1], fstBlast.seq2pos[1], scdBlast.seq2pos[0],seq2dtce, 'type']))

    #Sort the list by seq1pos1
    gapList.sort(key=lambda x : x[1])
    #create gap family, then add it to the master list
    if len(gapList) > 0:
        newGapFamily = GapFamily(gapList)
        masterFamilyList.append(newGapFamily)

def getRefFastaName(refFasta):
    return SeqIO.read(refFasta, 'fasta').name

def createCombinedFasta(fastaList):
    loggerMain.info('Combining fastas')
    with open('combinedFasta.temp', 'w') as output:
        for fasta in fastaList:
            for record in SeqIO.parse(fasta, 'fasta'):
                if (str(record.name), fasta) not in fileList:
                    fileList.append((str(record.name), fasta))
                if record.name not in fastaListLengths.keys():
                    fastaListLengths[record.name] = len(record.seq)
                output.write('>' + str(record.name) + '\n')
                output.write(str(record.seq) + '\n')
    loggerMain.info('Fastas combined')

def removeBlastDBFiles():
    loggerMain.info('cleaning Blast DB files...')
    dirFiles = os.listdir(os.curdir)
    removedFileList = []
    removedFileCounter = 0
    for file in dirFiles:
        if file.split('.')[0] == 'dbTemp':
            removedFileCounter += 1
            removedFileList.append(file)
            os.remove(file)
    loggerMain.info('{} files removed: {}'.format(removedFileCounter, removedFileList))

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
            findAllGaps(family, thresholdMin=minThres, thresholdAnalysis=maxThres, masterFamilyList=masterFamilyList)

def equalizeBorderSizes(gap1borders, gap2borders):

    # remember: getBorderSequence returns a tuple: ((leftBorderSize, leftBorderSeq), (rightBorderSize, rightBorderSeq))
    # borders might have different sizes, so we have to equalize the border size
    # start with the left border:
    if (gap1borders[0][0] + gap2borders[0][0]) > 0:
        if gap1borders[0][0] != gap2borders[0][0]:
            loggerMain.info('left side borders are not of the same size! gap1: {}, gap2: {}'.format(gap1borders[0][0], gap2borders[0][0]))

            if gap1borders[0][0] > gap2borders[0][0]:
                #we calculate the difference with the total bordar size because it has the same result
                diff = gap1borders[0][0] - gap2borders[0][0]
                gap1borders[0] = (gap1borders[0][0] - diff, gap1borders[0][1][diff:])
                loggerMain.info('Fixed gap1 border: now gap1size is {} and gap2size is {}'.format(gap1borders[0][0], gap2borders[0][0]))
            elif gap2borders[0][0] > gap1borders[0][0]:
                diff = gap2borders[0][0] - gap1borders[0][0]
                gap2borders[0] = (gap2borders[0][0] - diff, gap2borders[0][1][diff:])
                loggerMain.info(
                    'Fixed gap2 border: now gap1size is {} and gap2size is {}'.format(gap1borders[0][0], gap2borders[0][0]))
            loggerMain.debug('gap1leftBorder: {},  gap2leftBorder: {}'.format(len(gap1borders[0][1]), len(gap2borders[0][1])))


    #right border:
    if (gap1borders[1][0] + gap2borders[1][0]) > 0:
        if gap1borders[1][0] != gap2borders[1][0]:
            loggerMain.debug('right side borders are not of the same size! gap1: {}, gap2: {}'.format(gap1borders[1][0],
                                                                                                    gap2borders[1][0]))
            if gap1borders[1][0] > gap2borders[1][0]:
                diff = gap1borders[1][0] - gap2borders[1][0]
                gap1borders[1] = (gap1borders[1][0] - diff, gap1borders[1][1][:len(gap1borders[1][1]) - diff])
                loggerMain.info(
                    'Fixed gap1 border: now gap1size is {} and gap2size is {}'.format(gap1borders[1][0], gap2borders[1][0]))
            elif gap2borders[1][0] > gap1borders[1][0]:
                diff = gap2borders[1][0] - gap1borders[1][0]
                gap2borders[1] = (gap2borders[1][0] - diff, gap2borders[1][1][:len(gap2borders[1][1]) - diff])
                loggerMain.info(
                    'Fixed gap2 border: now gap1size is {} and gap2size is {}'.format(gap1borders[1][0], gap2borders[1][0]))
            loggerMain.debug('gap1rightBorder: {}, gap2rightBorder: {}'.format(len(gap2borders[1][1]), len(gap2borders[1][1])))

        else:
            loggerMain.info('No need to fix borders: gap1size is {} and gap2size is {}'.format(gap1borders[1][0], gap2borders[1][0]))

def compareTwoGaps(gap1, gap2, size, storeResults = True):
    # Get the border sequences for both gaps
    gap1borders = gap1.getBorderSequence(size)
    gap2borders = gap2.getBorderSequence(size)
    loggerMain.info(gap1borders)
    loggerMain.info(gap2borders)

    #equalize borders if required
    if gap1borders[0][1] != None and gap2borders[0][1] != None and gap1borders[0][0] != gap2borders[0][0]:
        equalizeBorderSizes(gap1borders, gap2borders)
    if gap1borders[1][1] != None and gap2borders[1][1] != None and gap1borders[1][0] != gap2borders[1][0]:
        equalizeBorderSizes(gap1borders, gap2borders)



    #call stretcher and get the alignment. first, prepare the fasta files
    tempFastaFiles = ['./temp/gap1left.temp.fasta', './temp/gap1right.temp.fasta', './temp/gap2left.temp.fasta', './temp/gap2right.temp.fasta']
    rightalnResults = None
    leftalnResults = None
    loggerMain.info('comparing left border')
    if gap1borders[0][1] != None and gap2borders[0][1] != None:
        with open(tempFastaFiles[0], 'w') as gap1left:
            gap1left.write('>gap1left\n')
            gap1left.write(str(gap1borders[0][1]))
        with open(tempFastaFiles[2], 'w') as gap2left:
            gap2left.write('>gap2left\n')
            gap2left.write(str(gap2borders[0][1]))
        leftalnResults = alignSeqsByBlast('./temp/gap1left.temp.fasta', './temp/gap2left.temp.fasta', gap1borders[0][0])
    else:
        loggerMain.info('No sequences for the left border: border 1 is {} and border 2 is {}'.format(len(gap1borders[0][1]), len(gap2borders[0][1])))
        rightalnResults = {'Identity': -1, 'Gaps': 0, 'Matches%': 0, 'AlnLen': 0, 'Mismatches': 0}

    loggerMain.info('comparing right border')
    if gap1borders[1][1] != None and gap2borders[1][1] != None:
        with open(tempFastaFiles[1], 'w') as gap1right:
            gap1right.write('>gap1right\n')
            gap1right.write(str(gap1borders[1][1]))
        with open(tempFastaFiles[3], 'w') as gap2right:
            gap2right.write('>gap2right\n')
            gap2right.write(str(gap2borders[1][1]))
        rightalnResults = alignSeqsByBlast('./temp/gap1right.temp.fasta', './temp/gap2right.temp.fasta', gap1borders[1][0])
    else:
        loggerMain.info('No sequences for the right border')
        rightalnResults = {'Identity': -1, 'Gaps': 0, 'Matches%': 0, 'AlnLen': 0, 'Mismatches': 0}


    #check if they are similar. If they are, then add borders
    response = [False, False]
    #todo - check if reverse complement sequences fuck with identity(use blast then)
    #All of this is code for stretcher
    '''
    if leftalnResults['Identity'] > 95 and leftalnResults['Gaps'] < 1:
        loggerMain.debug('left sides match')
        response[0] = True
        if storeResults == True:
            gap1.addBorderGapRef(gap2, 'left')
            gap2.addBorderGapRef(gap1, 'left')

    if rightalnResults['Identity'] > 95 and rightalnResults['Gaps'] < 1:
        loggerMain.debug('right sides match')
        response[0] = True
        if storeResults == True:
            gap1.addBorderGapRef(gap2, 'right')
            gap2.addBorderGapRef(gap1, 'right')
    '''
    loggerMain.debug('leftalnResults: {}'.format(leftalnResults))
    if leftalnResults['Matches%'] > 95 and leftalnResults['AlnLen'] > 0.95:
        loggerMain.debug('left sides match')
        response[0] = True
        if storeResults == True:
            gap1.addBorderGapRef(gap2, 'left')
            gap2.addBorderGapRef(gap1, 'left')

    loggerMain.debug('rightalnResults: {}'.format(rightalnResults))
    if rightalnResults['Matches%'] > 95 and leftalnResults['AlnLen'] > 0.95:
        loggerMain.debug('right sides match')
        response[1] = True
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
    totalNoOfCombinations = len(gapFamily1.getGapsByType('analyze')) * len(gapFamily2.getGapsByType('analyze'))
    currComb = 0
    combList = iter([(x, y) for x in gapFamily1.getGapsByType('analyze') for y in gapFamily2.getGapsByType('analyze')])

    if not os.path.exists('./temp/'):
        os.mkdir('temp')

    for gapPair in combList:
        currComb += 1
        print('{} / {}'.format(currComb, totalNoOfCombinations))
        if gapPair[0].family != gapPair[1].family:
            loggerMain.debug('comparing gap pair {} - {}'.format(gapPair[0].name, gapPair[1].name))
            compareTwoGaps(gapPair[0], gapPair[1], size)

def checkForGapClusters(masterFamilyList):
    gapClusterList = []
    gapMasterList = []
    gapCheckedList = []
    #populate gapMasterList
    for gapFamily in masterFamilyList:
        print(len(gapFamily.getGapsByType('analyze')))
        for gap in gapFamily.getGapsByType('analyze'):
            gapMasterList.append(gap)

    #find clusters
    for gap in gapMasterList:
        if gap in gapCheckedList:
            continue
        else:
            print('looking for partners for {}'.format(gap.name))
            if set(gap.leftMatch) == set(gap.rightMatch):
                newCluster = [gap]
                for borderGap in gap.leftMatch:
                    if set(borderGap.leftMatch) == set(borderGap.rightMatch):
                        print('\t1/2')
                        newList = copy.copy(borderGap.leftMatch)
                        newList.remove(gap)
                        newList.append(borderGap)
                        if set(gap.leftMatch) == set(newList):
                            print('\t2/2')
                            newCluster.append(borderGap)
                for gap in newCluster:
                    gapCheckedList.append(gap)
                gapClusterList.append(newCluster)
    return gapClusterList

def dumpClusterData(gapClusterList):
    # get a list of all sequences analyzed:
    analyzedList = []
    for filePair in fileList:
        if filePair[1] != refFasta:
            analyzedList.append(filePair[0])

    #also get a list for all gaps not part of a gapCluster
    unclusteredGaps = []
    # create and populate the fastaList
    fastaList = []
    for file in fileList:
        if file[1] != refFasta:
            fastaList.append(file[1])

    performProdigalTraining(fastaList)
    #If the cluster is a real cluster (it has more than one gap), then add it to the unclustered list
    for i, gapCluster in enumerate(gapClusterList):
        #We want to have in the info file which sequences do not have the
        if len(gapCluster) > 1:
            clusterGapNames = []
            clusterGapIndexes = []
            orfDict = getORFs(gapCluster, gapClusterList.index(gapCluster))

            for gap in gapCluster:
                clusterGapNames.append(gap.parents[0])
                clusterGapIndexes.append(gapCluster.index(gap))
            infoFile = open('./other/gapCluster_{}.txt'.format(i), 'w')
            infoFile.write('NAME\tPOS1\tPOS2\tLEN\tNOOFGENES\n')
            fastaFile = open('./other/gapCluster_{}.fasta'.format(i), 'w')
            for name in analyzedList:
                if name in clusterGapNames:
                    index = clusterGapNames.index(name)
                    gap = gapCluster[index]
                    # write fasta file:
                    fastaFile.write('>{}-{}\n'.format(gap.name, gap.parentPos[1] - gap.parentPos[0]))
                    fastaFile.write('{}\n'.format(gap.getParentSeq()[gap.parentPos[0]:gap.parentPos[1]]))
                    # write info file:
                    infoFile.write('{}\t{}\t{}\t{}\t{}\n'.format(gap.parents[0], gap.parentPos[0], gap.parentPos[1],
                                                             gap.parentPos[1] - gap.parentPos[0], len(orfDict[name])))
                else:
                    # write info file:
                    infoFile.write('{}\t{}\t{}\t{}\n'.format(name, 'NULL', 'NULL', 'NULL'))

            infoFile.close()
            fastaFile.close()
        else:
            unclusteredGaps.append(gapCluster[0])

def performProdigalTraining(list):
    createCombinedFasta(list)
    subprocess.call([prodigalPath+'prodigal', '-p', 'train', '-i', 'combinedFasta.temp', '-t', 'training.tmp.trn'])
    os.remove('combinedFasta.temp')
    return 'training.tmp.trn'

def getORFs(gapCluster, index):
    seqDict = {}
    for gap in gapCluster:
        name = 'gapCluster_{}_{}'.format(index, gap.parents[0])
        seq = gap.getParentSeq()[gap.parentPos[0]:gap.parentPos[1]]
        seqDict[name] = seq

    with open('combinedGapFastas.temp.fasta', 'w') as output:
        for key in seqDict:
            output.write('>{}\n'.format(key))
            output.write('{}\n'.format(seqDict[key]))

    subprocess.call([prodigalPath+'prodigal', '-i', 'combinedGapFastas.temp.fasta', '-p', 'single', '-t', 'training.tmp.trn', '-f', 'sco', '-o', 'outputGaps_{}.sco'.format(index)])
    orfDict = parseScoFile('outputGaps_{}.sco'.format(index))
    os.remove('outputGaps_{}.sco'.format(index))
    return orfDict

def parseScoFile(scoFile):
    scoDict = {}
    with open(scoFile, 'r') as scoInput:
        newKey = None
        while True:
            line = scoInput.readline()
            if not line:
                break
            if line[0:3] == '# S':
                newKey = line.split(';')[-1].split('"')[-2].split('_')[-1]
                scoDict[newKey] = []
            elif line[0] == '>':
                splitLine = line.split('_')
                geneInfo = (splitLine[1], splitLine[2], splitLine[3].rstrip('\n'))
                print(geneInfo)
                scoDict[newKey].append(geneInfo)
    return scoDict


if __name__ == '__main__':


    blastFamilies = performBlastAgainstCore(fastaList, refFasta, 500, 2000, masterFamilyList)
    os.remove(os.curdir + '/blastSeqs.blastn')
    for family in masterFamilyList:
        family.equalize()

    #Perform comparisons
    familyCombList = itertools.combinations(masterFamilyList, 2)
    for familyPair in familyCombList:
        setup_logger('herringFamily', 'herring.{}.{}.log'.format(familyPair[0].parents[0], familyPair[1].parents[0]))
        loggerMain = logging.getLogger('herringFamily')
        compareGapFamilies(familyPair[0], familyPair[1], 5000)

    #go back to old log here
    setup_logger('herringMain', 'herring.main.log')
    loggerMain = logging.getLogger('herringMain')
    #Save the gap file, just in case something breaks so we don't have to do the analysis all over again

    saveHerringData('masterListTest.pk1')
    loadedStuff = loadHerringData('masterListTest.pk1')
    masterFamilyList = loadedStuff[0]
    fileList = loadedStuff[1]


    for family in masterFamilyList:
        family.writeGapInfo(family.findOrphans(type = ('left', 'right')), 'diagnostic-All' + '-' + str(family.parents[0]) + '.txt')
        family.writeGapInfo([item for item in family.gapList if item not in family.findOrphans(type=('left', 'right', 'both'))], 'diagnostic-Good' + '-' + str(family.parents[0]) + '.txt')

    gapClusterList = checkForGapClusters(masterFamilyList)
    dumpClusterData(gapClusterList)






