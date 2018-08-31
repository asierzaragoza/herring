#todo - right now the fasta has to have the same name as its header to work. Not a critical error, but a pain in the ass

import subprocess, fnmatch, os, copy, logging, itertools, sys, pickle, time, math
from builtins import print

from Bio import SeqIO
import multiprocessing

#Log setup so we can have more than 1 log at the same time (shamelessly stolen from stackOverflow)
def setup_logger(logger_name, log_file, level=logging.DEBUG):
    l = logging.getLogger(logger_name)
    formatter = logging.Formatter('%(message)s')
    fileHandler = logging.FileHandler(log_file, mode='a')
    fileHandler.setFormatter(formatter)
    l.setLevel(level)
    l.addHandler(fileHandler)

setup_logger('herringMain', 'herring.main.log')
loggerMain = logging.getLogger('herringMain')

sys.path.insert(0, '/home/asier/PycharmProjects/flex2')
import blastParser

#Windows stretcher path
#stretcherPath = 'C:\mEMBOSS\stretcher.exe'
stretcherPath = '/usr/bin/stretcher'
#blastPath = '/opt/ncbi-blast-2.6.0+/bin/'
blastPath = '/home/rohit/alberto/ncbi-blast-2.7.1+/bin/'
prodigalPath = '/home/rohit/asier/programs/prodigal/'

###Lactobacillus stuff
#fastaList = ['B21.fasta', 'C410L1.fasta', 'DOMLa.fasta', 'JBE245.fasta', 'JBE490.fasta',
#             'LP2.fasta', 'LP3.fasta', 'LY78.fasta', 'LZ227.fasta', 'MF1298.fasta', 'RI113.fasta',
#             'SRCM100434.fasta', 'ZJ316.fasta', 'WCFS1.fasta']
#fastaList = ['WCFS1.fasta', 'LZ227.fasta', 'B21.fasta']
#refFasta = 'lactobacillus.coreWCFS1.fasta'

###Methanosarcina stuff
#fastaList = ['C16.fasta', 'Goe1.fasta', 'LYC.fasta', 'S-6.fasta', 'SarPi.fasta', 'Tuc01.fasta', 'WWM610.fasta']
fastaList = ['S-6.fasta', 'LYC.fasta', 'WWM610.fasta']
refFasta = 'methanosarcina.core.fasta'

###Clostridioides stuff
#fastaList = ['630.fasta', '630Drm.fasta', '630Dta.fasta', 'BI9.fasta', 'CD196.fasta', 'CD630DERM.fasta', 'CF5.fasta',
#             'M68.fasta', 'M120.fasta', 'NAP11.fasta', 'R0104a.fasta', 'R20291.fasta', 'W0003a.fasta',
#             'W0022a.fasta', 'W0023a.fasta']
#refFasta = 'clostridioides.core.fasta'

###Listeria stuff
#fastaList = ['CFSAN023459.fasta', 'EGD.fasta', 'EGDe.fasta', 'FDA00011238.fasta', 'FORC049.fasta', 'H34.fasta',
#             'J0161.fasta', 'L2624.fasta', 'L2676.fasta', 'LM11.fasta', 'LM850658.fasta', 'M7.fasta',
#             'NCTC10357.fasta', 'NTSN.fasta', 'WSLC1019.fasta']
#refFasta = 'listeria.core.fasta'


fastaListLengths = {}
fileList = []
masterFamilyList = []

#How many cores will be used?
cores = 3
#Hopefully on a near future user will enter this number, but today is not that day. Today you must write a number

#Stuff anti-weird numbers
if cores < 1:
    print('Not a good number of cores. Please insert a natural number')
    #Natural numbers are those positive ones excluding zero and without decimals. Dont change it, looks cool
    exit(-1)
f = math.factorial
max_number_of_combinations = f(cores) // f(2) // f(cores - 2)
#author's note: we use the integer division (//) to avoid possible oveflows
if cores > max_number_of_combinations:
    cores = max_number_of_combinations


class GapFamily:
    def __init__(self, gapList):
        self.parents = (gapList[0][0], gapList[0][4])
        self.gapList = []
        for gapItem in gapList:
            newGap = Gap(gapItem, self)
            self.gapList.append(newGap)

    def writeGapInfo(self, gapList, filename):
        #Super duper ultra cool new mode to check and create (if needed) the directory for the filename
        # Proudly made in StackOverflow
        os.makedirs(os.path.dirname('./logs/' + filename), exist_ok=True)
        with open('./logs/' + filename, "w") as output:
            output.write('NUMBER OF GAPS: {}\n\n'.format(len(gapList)))
            for i, gap in enumerate(gapList):
                gapStr = gap.buildInfoStr()
                output.write(gapStr)
                if i < len(gapList) - 1:
                    output.write('-----------------------------------\n\n')

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

    def findAdjacentGaps(self, gap):
        #sort gaps by seq1pos
        self.gapList.sort(key=lambda gap: gap.parentPos[0])
        #get the index for the target gap
        targetIndex = self.gapList.index(gap)
        prevGap = None
        nextGap = None
        loggerMain.info('target index: {}, gapListLen: {}'.format(targetIndex, len(self.gapList)))
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
            loggerMain.info('fixed: parent[0] (parent) is {} and parent[1] (core) is {}'.format(self.parents[0],
                                                                                                self.parents[0]))
        for gap in self.gapList:
            loggerMain.info('gap.parents[0] = {}, family.parents[0] is {}'.format(gap.parents[0], self.parents[0]))
            if gap.parents[0] != self.parents[0]:

                newParent = copy.copy(gap.corePos)
                newCore = copy.copy(gap.parentPos)

                gap.parents = self.parents
                gap.corePos = newCore
                gap.parentPos = newParent
            loggerMain.info('parents = {} , corePos = {} , parentPos = {}'.format(gap.parents, gap.corePos,
                                                                                  gap.parentPos))

    def getGapsByType(self, type='analyze'):
        targetList = []
        for gap in self.gapList:
            if gap.type == type:
                targetList.append(gap)
        return targetList

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
        self.gapLength = self.corePos[1] - self.corePos[0]

        if self.gapLength < 0:
            self.gapLength *= -1
            self.corePos = (self.corePos[1], self.corePos[0])

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
        #Check if we can pick a sequence of the required size without taking a gap by mistake:
        adjGaps = self.family.findAdjacentGaps(self)
        #Let's start with the left border
        if adjGaps[0] is not None:
            leftGapPos = adjGaps[0].parentPos[1]
            if (self.parentPos[0] - leftGapPos) < size:
                loggerMain.info('there is not enough space for a full left border! size: {}, space between gaps: {}'
                                 ''.format(size, self.parentPos[0] - leftGapPos))
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
                loggerMain.info('there is not enough space for a full right border! size: {}, space between gaps: {}'
                                 ''.format(size,rightGapPos - self.parentPos[1]))
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
        loggerMain.info('borders for {} found: left border is {}, while right border is {}'.format(self.name,
                                                        [(self.parentPos[0] - self.leftBorder[0]), self.parentPos[0]],
                                                        [self.parentPos[1], self.parentPos[1] + self.rightBorder[0]]))
        return [gapSeqLeft, gapSeqRight]

    def addBorderGapRef(self, gap, type):
        if type == 'left':
            self.leftMatch.append(gap)
        elif type == 'right':
            self.rightMatch.append(gap)

    def getOrphanList(self):
        orphanList = []
        for leftGap in self.leftMatch:
            if leftGap not in self.rightMatch:
                orphanList.append((leftGap, 'left'))

        for rightGap in self.rightMatch:
            if rightGap not in self.leftMatch:
                orphanList.append((rightGap, 'right'))
        return orphanList

    def isGapOrphan(self):
        if len(self.leftMatch) == len(self.rightMatch) and len(self.leftMatch) + len(self.rightMatch) != 0:
            return [True, 'both']
        elif len(self.leftMatch) < len(self.rightMatch):
            return [True, 'right']
        else:
            return [True, 'left']
            #Notice that there are some gaps that will be classified as left even they are not, like those with len
            # for both borders equal zero, but we do this on purpose because the only call for this function is to
            # classify the gaps in diagnostic Good or All. Actually we could dessign this function only with the two
            # first lines and a return [True, 'left'] at the end, user wont notice any difference.
            # I wont swap the function as told just because other programmers may modify this section in order to expand
            # its use, then maybe they find useful this 'extended version'


def saveHerringData(filename):
    with open(filename, 'wb') as output:
        pickle.dump(masterFamilyList, output, pickle.HIGHEST_PROTOCOL)
        pickle.dump(fileList, output, pickle.HIGHEST_PROTOCOL)

def loadHerringData(filename):
    with open(filename, 'rb') as input:
        familylist = pickle.load(input)
        filelist = pickle.load(input)
    return(familylist, filelist)

def alignSeqsByBlast(seq1, seq2, borderSize, name):
    loggerMain.info('Aligning 2 Sequences by Blast')
    outfile = name + 'blastResults.temp.blastn'
    loggerMain.info('Running blastn...')
    subprocess.call([blastPath + 'blastn', '-query', seq1, '-subject', seq2, '-out', outfile, '-outfmt', '6'])
    statsDict = {'Matches%': 0, 'AlnLen': 0, 'Gaps': 0, 'Mismatches': 0}
    if os.stat(outfile).st_size > 0:
        loggerMain.info('There was a match!')
        with open(outfile, 'r') as blastHandle:
            statsDict = {'Matches%':0 , 'AlnLen': 0, 'Gaps': 0, 'Mismatches' : 0}
            #remember: queryID / subjectID / percentage of identical matches / alignment length / n of mismatches
            # / n of gaps / start in query / end in query / start in subject / end in subject / E-value / bitscore
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
        loggerMain.info('Iteration n {}'.format(i))
        # Get the 2 blasts to compare
        fstBlast = copy.copy(blastList[i])
        scdBlast = copy.copy(blastList[i + 1])
        # Check if the blasts are reversed, and if they are fix them
        if fstBlast.seq1pos[1] < fstBlast.seq1pos[0]:
            fstBlast.seq1pos = (fstBlast.seq1pos[1], fstBlast.seq1pos[0])
        if scdBlast.seq1pos[1] < scdBlast.seq1pos[0]:
            scdBlast.seq1pos = (scdBlast.seq1pos[1], scdBlast.seq1pos[0])
        # calculate distance between both blasts in both sequences
        loggerMain.info('blast1pos: {} - {}'.format(fstBlast.seq1pos[0], fstBlast.seq1pos[1]))
        loggerMain.info('blast2pos: {} - {}'.format(scdBlast.seq1pos[0], scdBlast.seq1pos[1]))
        seq1dtce = abs(scdBlast.seq1pos[0] - fstBlast.seq1pos[1])
        seq2dtce = abs(scdBlast.seq2pos[0] - fstBlast.seq2pos[1])
        loggerMain.info('Distances between blasts: seq1 {}, seq2 {}'.format(seq1dtce, seq2dtce))
        total = 0
        # Check if the distances in both sequences meet the required threshold
        if seq1dtce > thresholdMin:
            total += 1
        if seq2dtce > thresholdMin:
            total += 1
        # then assign a category depending on the distances
        # gapList is: parent1, seq1pos1, seq1pos2, seq1dtce, parent2, seq2pos1, seq2pos2, seq2dtce, type
        if total == 0:
            loggerMain.info('Total is {}, which means there is no gap'.format(total))
            pass
        elif total > 0:
            loggerMain.info('Total is {}, which means there is a gap'.format(total))
            type = 'ignore'

            if seq1dtce + seq2dtce > thresholdAnalysis:
                type = 'analyze'

            gapList.append([fstBlast.parents[0], fstBlast.seq1pos[1], scdBlast.seq1pos[0], seq1dtce,
                            fstBlast.parents[1], fstBlast.seq2pos[1], scdBlast.seq2pos[0], seq2dtce, type])
            loggerMain.info('Sotring gap as gapList {}'.format(
                [fstBlast.parents[0], fstBlast.seq1pos[1], scdBlast.seq1pos[0], seq1dtce, fstBlast.parents[1],
                 fstBlast.seq2pos[1], scdBlast.seq2pos[0], seq2dtce, type]))

    # Sort the list by seq1pos1
    gapList.sort(key=lambda x: x[1])
    # create gap family, then add it to the master list
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

        blastHits = blastParser.parseBlastFile('blastSeqs.blastn', minIdentity=85, minAln=1500)
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
            loggerMain.info('left side borders are not of the same size! gap1: {}, gap2: {}'.format(gap1borders[0][0],
                                                                                                    gap2borders[0][0]))
            if gap1borders[0][0] > gap2borders[0][0]:
                #we calculate the difference with the total bordar size because it has the same result
                diff = gap1borders[0][0] - gap2borders[0][0]
                gap1borders[0] = (gap1borders[0][0] - diff, gap1borders[0][1][diff:])
                loggerMain.info('Fixed gap1 border: now gap1size is {} and gap2size is {}'.format(gap1borders[0][0],
                                                                                                  gap2borders[0][0]))
            elif gap2borders[0][0] > gap1borders[0][0]:
                diff = gap2borders[0][0] - gap1borders[0][0]
                gap2borders[0] = (gap2borders[0][0] - diff, gap2borders[0][1][diff:])
                loggerMain.info('Fixed gap2 border: now gap1size is {} and gap2size is {}'
                                ''.format(gap1borders[0][0], gap2borders[0][0]))
            loggerMain.info('gap1leftBorder: {},  gap2leftBorder: {}'.format(len(gap1borders[0][1]),
                                                                              len(gap2borders[0][1])))
    #right border:
    if (gap1borders[1][0] + gap2borders[1][0]) > 0:
        if gap1borders[1][0] != gap2borders[1][0]:
            loggerMain.info('right side borders are not of the same size! gap1: {}, gap2: {}'.format(gap1borders[1][0],
                                                                                                    gap2borders[1][0]))
            if gap1borders[1][0] > gap2borders[1][0]:
                diff = gap1borders[1][0] - gap2borders[1][0]
                gap1borders[1] = (gap1borders[1][0] - diff, gap1borders[1][1][:len(gap1borders[1][1]) - diff])
                loggerMain.info(
                    'Fixed gap1 border: now gap1size is {} and gap2size is {}'.format(gap1borders[1][0],
                                                                                      gap2borders[1][0]))
            elif gap2borders[1][0] > gap1borders[1][0]:
                diff = gap2borders[1][0] - gap1borders[1][0]
                gap2borders[1] = (gap2borders[1][0] - diff, gap2borders[1][1][:len(gap2borders[1][1]) - diff])
                loggerMain.info('Fixed gap2 border: now gap1size is {} and gap2size is {}'.format(gap1borders[1][0],
                                                                                                  gap2borders[1][0]))
            loggerMain.info('gap1rightBorder: {}, gap2rightBorder: {}'.format(len(gap2borders[1][1]),
                                                                               len(gap2borders[1][1])))
        else:
            loggerMain.info('No need to fix borders: gap1size is {} and gap2size is {}'.format(gap1borders[1][0],
                                                                                               gap2borders[1][0]))

def compareTwoGaps(gap1, gap2, size, actualLog, name, storeResults = True):
    # Get the border sequences for both gaps
    gap1borders = gap1.getBorderSequence(size)
    gap2borders = gap2.getBorderSequence(size)
    actualLog.info(gap1borders)
    actualLog.info(gap2borders)
    #equalize borders if required
    if gap1borders[0][1] != None and gap2borders[0][1] != None and gap1borders[0][0] != gap2borders[0][0]:
        equalizeBorderSizes(gap1borders, gap2borders)
    if gap1borders[1][1] != None and gap2borders[1][1] != None and gap1borders[1][0] != gap2borders[1][0]:
        equalizeBorderSizes(gap1borders, gap2borders)
    #call stretcher and get the alignment. first, prepare the fasta files
    tempFastaFiles = [name + 'gap1left.temp.fasta', name + 'gap1right.temp.fasta',
                          name + 'gap2left.temp.fasta', name + 'gap2right.temp.fasta']

    actualLog.info('comparing left border')
    if gap1borders[0][1] != None and gap2borders[0][1] != None:
        with open(tempFastaFiles[0], 'w') as gap1left:
            gap1left.write('>gap1left\n')
            gap1left.write(str(gap1borders[0][1]))
        with open(tempFastaFiles[2], 'w') as gap2left:
            gap2left.write('>gap2left\n')
            gap2left.write(str(gap2borders[0][1]))
        leftalnResults = alignSeqsByBlast(name + 'gap1left.temp.fasta', name + 'gap2left.temp.fasta',
                                          gap1borders[0][0], name)
    else:
        actualLog.info('No sequences for the left border: border 1 is {} and border 2 is {}'
                       ''.format(len(gap1borders[0][1]), len(gap2borders[0][1])))
        leftalnResults = {'Identity': -1, 'Gaps': 0, 'Matches%': 0, 'AlnLen': 0, 'Mismatches': 0}

    actualLog.info('comparing right border')
    if gap1borders[1][1] != None and gap2borders[1][1] != None:
        with open(tempFastaFiles[1], 'w') as gap1right:
            gap1right.write('>gap1right\n')
            gap1right.write(str(gap1borders[1][1]))
        with open(tempFastaFiles[3], 'w') as gap2right:
            gap2right.write('>gap2right\n')
            gap2right.write(str(gap2borders[1][1]))
        rightalnResults = alignSeqsByBlast(name + 'gap1right.temp.fasta', name + 'gap2right.temp.fasta',
                                           gap1borders[1][0], name)
    else:
        actualLog.info('No sequences for the right border')
        rightalnResults = {'Identity': -1, 'Gaps': 0, 'Matches%': 0, 'AlnLen': 0, 'Mismatches': 0}

    #remove all temp fasta files
    for file in tempFastaFiles:
        if file in os.listdir(os.curdir):
            os.remove(file)
    #check if they are similar. If they are, then add borders
    response = [False, False]
    data = []
    actualLog.info('leftalnResults: {}'.format(leftalnResults))
    if leftalnResults['Matches%'] > 89 and leftalnResults['AlnLen'] > 0.89:
        actualLog.info('left sides match')
        response[0] = True
        if storeResults == True:
            gap1.addBorderGapRef(gap2, 'left')
            gap2.addBorderGapRef(gap1, 'left')
            data.append('left')

    actualLog.info('rightalnResults: {}'.format(rightalnResults))
    if rightalnResults['Matches%'] > 89 and leftalnResults['AlnLen'] > 0.89:
        actualLog.info('right sides match')
        response[1] = True
        if storeResults == True:
            gap1.addBorderGapRef(gap2, 'right')
            gap2.addBorderGapRef(gap1, 'right')
            data.append('right')

    if response[0] != response[1]:
        loggerMain.info('both borders do not match!, left border: {}, right border: {}'.format(response[0],
                                                                                                response[1]))
        if closeGapByCorePos(gap1, gap2, size, name):
            data.append('right')
    return (gap1, gap2, data)

def compareGapFamilies(gapFamily1, gapFamily2, size, actualLog):
    # we want all combinations, not permutations (permutations also include the order of the elements - so AB & BA)
    # combList = itertools.combinations(gapFamily1.gapList + gapFamily2.gapList, 2)

    # Actually we don't want all combinations either, we just want all combinations between 2 sets (e.g. if list1 = 25
    # and list2 = 15, there would be 780 combinations but only 300 combinations between both sets -> Try list
    # comprehensions instead, and build an iterator if it becomes too much of a memory burden (it will)
    currComb = 0
    combList = iter([(x, y) for x in gapFamily1.getGapsByType('analyze') for y in gapFamily2.getGapsByType('analyze')])

    familyName = './' + gapFamily1.parents[0] + ',' + gapFamily2.parents[0]
    if not os.path.exists('./temp/'):
        os.mkdir('./temp/')

    result = []
    for gapPair in combList:
        currComb += 1
        if gapPair[0].family != gapPair[1].family:
            actualLog.info('comparing gap pair {} - {}'.format(gapPair[0].name, gapPair[1].name))
            result.append(compareTwoGaps(gapPair[0], gapPair[1], size, actualLog, name = familyName))

    return result

def checkForGapClusters(masterFamilyList):
    gapClusterList = []
    gapMasterList = []
    gapCheckedList = []
    #populate gapMasterList
    for gapFamily in masterFamilyList:
        for gap in gapFamily.getGapsByType('analyze'):
            gapMasterList.append(gap)
    #find clusters
    for gap in gapMasterList:
        if gap not in gapCheckedList:
            if set(gap.leftMatch) == set(gap.rightMatch):
                newCluster = [gap]
                for borderGap in gap.leftMatch:
                    if set(borderGap.leftMatch) == set(borderGap.rightMatch):
                        newList = copy.copy(borderGap.leftMatch)
                        newList.remove(gap)
                        newList.append(borderGap)
                        if set(gap.leftMatch) == set(newList):
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
            orfDict = getORFs(gapCluster, i)

            for gap in gapCluster:
                clusterGapNames.append(gap.parents[0])
                clusterGapIndexes.append(gapCluster.index(gap))
            infoFile = open('./other/gapCluster_{}.txt'.format(i), 'w')
            infoFile.write('NAME\tPOS1\tPOS2\tLEN\tNOOFGENES\n')
            fastaFile = open('./other/gapCluster_{}.fasta'.format(i), 'w')
            for name in analyzedList:
                if name in orfDict.keys():
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

def closeGapByCorePos(gap1, gap2, size, name, storeResults = True):
    if (gap1.corePos[1] - gap2.corePos[1]) > 0:
        longGap = gap1
        targetGap = gap2
        coreDiff = gap1.corePos[1] - gap2.corePos[1]
    else:
        longGap = gap2
        targetGap = gap1
        coreDiff = gap2.corePos[1] - gap1.corePos[1]
    longBorder = longGap.getBorderSequence(size + coreDiff)[1]
    targetBorder = targetGap.getBorderSequence(size + coreDiff)[1]
    if longBorder[0] < coreDiff:
        loggerMain.info('\tNot enough space for the check! long border size: {}, coreDiff: {}'.format(longBorder[0],
                                                                                                      coreDiff))
        return False
    else:
        targetBorder = [targetBorder[0] - coreDiff, targetBorder[1][coreDiff:]]
        #longBorder = [longBorder[0] - coreDiff, longBorder[1][coreDiff:]]
        eqBorders = equalizeBorders(longBorder, targetBorder)
        tempFastaFiles = [name + 'border1Reg.temp.fasta', name + 'border2Reg.temp.fasta']

        loggerMain.info('comparing right border')
        if eqBorders[0][1] != None and eqBorders[1][1] != None:
            with open(tempFastaFiles[0], 'w') as border1fasta:
                border1fasta.write('>border1\n')
                border1fasta.write(str(eqBorders[0][1]))
            with open(tempFastaFiles[1], 'w') as border2fasta:
                border2fasta.write('>border2\n')
                border2fasta.write(str(eqBorders[1][1]))
            blastResults = alignSeqsByBlast(tempFastaFiles[0], tempFastaFiles[1], eqBorders[0][0], name)
        else:
            loggerMain.info('No sequences for the left border: border 1 is {} and border 2 is {}'.format(
                len(eqBorders[0][1]), len(eqBorders[0][1])))
            blastResults = {'Identity': -1, 'Gaps': 0, 'Matches%': 0, 'AlnLen': 0, 'Mismatches': 0}

        loggerMain.info('\t{}'.format(blastResults))
        if blastResults['Matches%'] > 80 and blastResults['AlnLen'] > 0.90:
            if storeResults == True:
                gap1.addBorderGapRef(gap2, 'right')
                gap2.addBorderGapRef(gap1, 'right')
                return True
    return False

def equalizeBorders(border1, border2, typeOfBorder='right'):
    loggerMain.info('running equalizeBorders')
    if (border1[0] + border2[0]) > 0 and border1[0] == border2[0]:
        loggerMain.info('No need to equalize, both sequences are the same size. border 1 {}, border 2 {}'
                        ''.format(border1[0], border2[0]))
        return (border1, border2)
    elif (border1[0] + border2[0]) > 0 and typeOfBorder == 'right':
        loggerMain.info('sequences are not the same size, treating them as RIGHT sequences. border 1 {}, border 2 {}'
                        ''.format(border1[0], border2[0]))

        if border1[0] > border2[0]:
            diff = border1[0] - border2[0]
            border1 = (border1[0] - diff, border1[1][:len(border1[1]) - diff])
            loggerMain.info('Fixed border1: now border1 is {} and selfSize is {}'.format(border1[0],border2[0]))
            return (border1, border2)

        elif border2[0] > border1[0]:
            diff = border2[0] - border1[0]
            border2 = (border2[0] - diff, border2[1][:len(border2[1]) - diff])
            loggerMain.info('Fixed border2: now spacerSize is {} and selfSize is {}'.format(border1[0], border2[0]))
            return (border1, border2)
    else:
        loggerMain.info('Both sequences have no length')

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

    subprocess.call([prodigalPath+'prodigal', '-i', 'combinedGapFastas.temp.fasta', '-p', 'single', '-t',
                     'training.tmp.trn', '-f', 'sco', '-o', 'outputGaps_{}.sco'.format(index)])
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
                scoDict[newKey].append(geneInfo)
    return scoDict

def checkForProblemType2(masterFamilyList, actualLogger):
    size = 5000
    toRemoveGaps = []
    newGaps = []
    for family in masterFamilyList:
        gapsToRemoveList = []
        newGapsFamilyList = []
        for gap in family.gapList:
            actualLogger.info('comparing left border')
            orphanList = [x[0] for x in gap.getOrphanList() if x[1] == 'left']
            closedGapList = []
            for leftGap in orphanList:
                if (gap.gapLength) > (leftGap.gapLength):
                    longGap = gap
                    longDtce = gap.gapLength
                    targetGap = leftGap
                    targetDtce = leftGap.gapLength
                else:
                    longGap = leftGap
                    longDtce = leftGap.gapLength
                    targetGap = gap
                    targetDtce = gap.gapLength
                familyName = './' + gap.family.parents[0] + ',' + leftGap.family.parents[0]
                results = iterativeCloseGaps(targetGap, longGap, targetDtce, longDtce, size, actualLogger,
                                         name = familyName)
                if results is not None:
                    closedGapList.append(results)
            for res in closedGapList:
                gapsToRemoveList += res[0]
                newGapsFamilyList.append(res[1])

        toRemoveGaps += gapsToRemoveList
        for gap1 in newGapsFamilyList:
            equal = False
            for gap2 in newGaps:
                if set(gap1.parentPos) == set(gap2.parentPos) and set(gap1.corePos) == set(gap2.corePos) and \
                        set(gap1.family) == set(gap2.family):
                    equal = True
                    break
            if equal == False:
                newGaps.append(gap1)

    toRemoveGaps.sort(key=lambda x: x.name)
    cleanGapsToRemove = []
    for i in range(0, len(toRemoveGaps) - 1):
        if i < len(toRemoveGaps) - 2:
            results = toRemoveGaps[i].name == toRemoveGaps[i + 1].name
            if results is False:
                cleanGapsToRemove.append(toRemoveGaps[i])
        else:
            if toRemoveGaps[i] not in cleanGapsToRemove:
                cleanGapsToRemove.append(toRemoveGaps[i])
            if toRemoveGaps[i + 1] not in cleanGapsToRemove:
                cleanGapsToRemove.append(toRemoveGaps[i + 1])

    newGaps.sort(key=lambda x: x.name)
    actualLogger.info('purging gaps: gaps in cleanGapsToRemove: {}'.format(len(cleanGapsToRemove)))
    for family in masterFamilyList:
        actualLogger.info('checking family {}'.format(family.parents))
        removeList = []
        for i, gap in enumerate(family.gapList):
            if gap in cleanGapsToRemove:
                removeList.append(gap)
                actualLogger.info('\t{}\tgap {} marked for removal'.format(i, gap.name))
            else:
                actualLogger.info('\t{}\tgap not in removalList, check matchLists for mentions'.format(i))
                removalLeft = []
                for matchGap in gap.leftMatch:
                    if matchGap in cleanGapsToRemove:
                        removalLeft.append(matchGap)
                for markedGap in removalLeft:
                    gap.leftMatch.remove(markedGap)

                removalRight = []
                for matchGap in gap.rightMatch:
                    if matchGap in cleanGapsToRemove:
                        removalRight.append(matchGap)

                for markedGap in removalRight:
                    gap.rightMatch.remove(markedGap)
        for markedGap in removeList:
            family.gapList.remove(markedGap)

    gapsByName = []
    for gap in newGaps:
        gap.family.gapList.append(gap)
        gapsByName.append(gap.name)

    # Compare the new gaps to get new relations
    actualLogger.info('\n\nASSINGING BORDERS TO THE NEW GAPS\n\n')
    dangerousCombinations = []
    for family in masterFamilyList:
        combList = iter([(x, y) for x in newGaps for y in family.getGapsByType('analyze')])
        if not os.path.exists('./temp/'):
            os.mkdir('./temp/')

        for gapPair in combList:
            #There is the case in which two gaps are written from different families if both of them are in newSeqs
            # In order to prevent duplicity is the second and third condition. We'll let only one of those combinations
            # happen and the others will be discarted.
            if gapPair[0].family != gapPair[1].family and (gapPair[0].name, gapPair[1].name) not in \
                    dangerousCombinations and (gapPair[1].name, gapPair[0].name) not in dangerousCombinations:
                if gapPair[0].name in gapsByName and gapPair[1].name in gapsByName:
                    dangerousCombinations.append((gapPair[0].name, gapPair[1].name))
                actualLogger.info('comparing gap pair {} - {}'.format(gapPair[0].name, gapPair[1].name))
                familyName = './' + gapPair[0].family.parents[0] + ',' + gapPair[1].family.parents[0]
                compareTwoGaps(gapPair[0], gapPair[1], size, actualLogger, name=familyName)

def iterativeCloseGaps(targetGap, longGap, targetDtce, longDtce, size, actualLogger, name):
    targetGapList = [targetGap]
    iterativeTargetGap = targetGap
    iterativeDtce = targetDtce
    gapThreshold = 1000
    longDtce += gapThreshold

    while True:
        actualLogger.info('\tstarting iteration: iterativeDtce is {}, longDtce is {}, gapThreshold is {}'.format(
                                                                                iterativeDtce,longDtce, gapThreshold))
        adjIterativeTargetGap = iterativeTargetGap.family.findAdjacentGaps(iterativeTargetGap)[1]
        if adjIterativeTargetGap is None:
            actualLogger.info('\t No more gaps to the right, breaking loop')
            break
        actualLogger.info('\tadjIterativeTargetGap is {}, with parentPos of {}'.format(adjIterativeTargetGap.name,
                                                                           adjIterativeTargetGap.parentPos))

        actualLogger.info('\titerative dtce is {} + {}, longDtce + gapThreshold is ({} + {})'.format(iterativeDtce,
                                        adjIterativeTargetGap.gapLength, longDtce - gapThreshold, gapThreshold))
        iterativeDtce += adjIterativeTargetGap.gapLength

        if iterativeDtce > longDtce:
            actualLogger.info('\t iterative dtce > longDtce, try merging method')
            regResult = regularCloseGaps(targetGap, longGap, targetDtce, longDtce, size, actualLogger, name)
            if regResult is not None:
                return regResult
            else:
                actualLogger.info('regular method did not provide any results')
                return None

        # Get borders
        longBorder = longGap.getBorderSequence(size)[1]
        adjBorder = adjIterativeTargetGap.getBorderSequence(size)[1]

        eqBorders = equalizeBorders(longBorder, adjBorder)
        actualLogger.info('\teqBorders len: {}, {}'.format(len(eqBorders[0][1]), len(eqBorders[1][1])))

        # then run blast:
        tempFastaFiles = [name + 'border1Adj.temp.fasta', name + 'border2Adj.temp.fasta']

        actualLogger.info('comparing right border')
        if eqBorders[0][1] != None and eqBorders[1][1] != None:
            with open(tempFastaFiles[0], 'w') as border1fasta:
                border1fasta.write('>border1\n')
                border1fasta.write(str(eqBorders[0][1]))
            with open(tempFastaFiles[1], 'w') as border2fasta:
                border2fasta.write('>gap2left\n')
                border2fasta.write(str(eqBorders[1][1]))
            blastResults = alignSeqsByBlast(tempFastaFiles[0], tempFastaFiles[1], eqBorders[0][0], name)
        else:
            actualLogger.info('No sequences for the left border: border 1 is {} and border 2 is {}'.format(
                len(eqBorders[0][1]), len(eqBorders[0][1])))
            blastResults = {'Identity': -1, 'Gaps': 0, 'Matches%': 0, 'AlnLen': 0, 'Mismatches': 0}

        # analyzeBlastResults
        if blastResults['Matches%'] > 80 and blastResults['AlnLen'] > 0.90:
            # create a new gap object and return it + the gaps used to create it
            newGap = Gap([targetGap.parents[0], targetGap.parentPos[0], adjIterativeTargetGap.parentPos[1], 'None',
                 targetGap.parents[1], targetGap.corePos[0], adjIterativeTargetGap.corePos[1], 'None',
                 'analyze'], targetGap.family)

            newNo = int(round((float(targetGap.name.split('_')[-1]) +
                               float(adjIterativeTargetGap.name.split('_')[-1])) / 2, 0))
            newGap.name = str(str(newGap.parents[0]) + '-' + str(newGap.parents[1]) + '_' + str(newNo))
            targetGapList.append(adjIterativeTargetGap)
            return (targetGapList, newGap)

        else:
            targetGapList.append(adjIterativeTargetGap)
            iterativeTargetGap = adjIterativeTargetGap
            actualLogger.info('\t new iterativeGap is {}'.format(adjIterativeTargetGap.name))

def regularCloseGaps(targetGap, longGap, targetDtce, longDtce, size, actualLogger, name):
    actualLogger.info('\tAttempting to close gap via merging')
    longBorder = longGap.getBorderSequence(size)[1]
    actualLogger.info('\tlong dtce:{}, targetDtce:{}, size:{}'.format(longDtce, targetDtce, size))
    targetBorder = targetGap.getBorderSequence(size + (longDtce - targetDtce))[1]
    actualLogger.info('\tlong border: {}\n\ttargetBorder: {}, diffDtce: {}'.format(longBorder, targetBorder,
                                                                                   (longDtce - targetDtce)))
    if targetBorder[0] < (longDtce - targetDtce):
        actualLogger.info('\tNot enough space for the check')
        return None

    else:
        targetBorder = [targetBorder[0] - (longDtce - targetDtce), targetBorder[1][(longDtce - targetDtce):]]
        eqBorders = equalizeBorders(longBorder, targetBorder)
        # then run blast:
        tempFastaFiles = [name + 'border1Reg.temp.fasta', name + 'border2Reg.temp.fasta']

        actualLogger.info('comparing right border')
        if eqBorders[0][1] != None and eqBorders[1][1] != None:
            with open(tempFastaFiles[0], 'w') as border1fasta:
                border1fasta.write('>border1\n')
                border1fasta.write(str(eqBorders[0][1]))
            with open(tempFastaFiles[1], 'w') as border2fasta:
                border2fasta.write('>gap2left\n')
                border2fasta.write(str(eqBorders[1][1]))
            blastResults = alignSeqsByBlast(tempFastaFiles[0], tempFastaFiles[1], eqBorders[0][0], name)
        else:
            actualLogger.info('No sequences for the left border: border 1 is {} and border 2 is {}'.format(
                len(eqBorders[0][1]), len(eqBorders[0][1])))
            blastResults = {'Identity': -1, 'Gaps': 0, 'Matches%': 0, 'AlnLen': 0, 'Mismatches': 0}

        # analyzeBlastResults
        if blastResults['Matches%'] > 80 and blastResults['AlnLen'] > 0.90:
            newGap = Gap([targetGap.parents[0], targetGap.parentPos[0], targetGap.parentPos[1] + (longDtce -
                    targetDtce), 'None', targetGap.parents[1], targetGap.corePos[0], targetGap.corePos[1] + (longDtce
                    - targetDtce), 'None', 'analyze'], targetGap.family)
            newGap.name = newGap.name.split('_')[0] + '_' + targetGap.name.split('_')[-1] + 'M'
            return ([targetGap], newGap)
        else:
            return None

def worker(tuple):
    family1 = tuple[0]
    family2 = tuple[1]
    size = 5000
    #Create one log for each combination and pass the log's reference at the function's call
    setup_logger('herring_process', 'herring.{}.{}.log'.format(family1.parents[0], family2.parents[0]))
    actualLog = logging.getLogger('herring_process')
    actualLog.info("Now working on families {} and {}".format(family1.parents[0], family2.parents[0]))
    data = compareGapFamilies(family1, family2, size, actualLog)
    #data is a list of tuples, all combinations of gaps from the two families of the thread
    return data


if __name__ == '__main__':

    timeStart = time.time()
    #Magic number powered by last author
    thresholdAnalysis = 2000
    blastFamilies = performBlastAgainstCore(fastaList, refFasta, 500, thresholdAnalysis, masterFamilyList)
    os.remove(os.curdir + '/blastSeqs.blastn')
    for family in masterFamilyList:
        family.equalize()

    #go back to old log here
    setup_logger('herringMain', 'herring.main.log')
    loggerMain = logging.getLogger('herringMain')
    #Save the gap file, just in case something breaks so we don't have to do the analysis all over again

    #Perform comparisons
    familyCombList = itertools.combinations(masterFamilyList, 2)

    #initialize the thread pool
    loggerMain.info('Now launching the threads. Number of threads desired by user: {}'.format(cores))
    p=multiprocessing.Pool(cores)
    data = p.map(worker, familyCombList)
    #Wait until the threads have finished their work
    p.close()
    p.join()
    loggerMain.info('Threads have already finished (all of them, program is supposed to wait for them)')

    loggerMain.info('Now assigning each tuple element from the threads to corresponding family')
    for familyList in data:
        #Each element is a list of tuples which represents all the gaps between the two families
        for gapFound in familyList:
            #Each tuple is the data for the gaps. Its structure is (gap1, gap2, border)
            #Border is a list of strings. It can be empty, figure 'left', 'right' or both (but not the word 'both')
            gap1 = gapFound[0]
            gap2 = gapFound[1]
            border = gapFound[2]
            family1Name = gap1.family.parents
            family2Name = gap2.family.parents
            classFamily1 = None
            classFamily2 = None

            #Look for the gap's family
            for familyRoot in masterFamilyList:
                if classFamily1 == None and familyRoot.parents == family1Name:
                    classFamily1 = familyRoot
                if classFamily2 == None and familyRoot.parents == family2Name:
                    classFamily2 = familyRoot
                if classFamily1 != None and classFamily2 != None:
                    break

            #Now look for the gap whose name is the same as the gap's thread
            centralGap1 = None
            centralGap2 = None
            for candidateGap in classFamily1.gapList:
                if candidateGap.name == gap1.name:
                    centralGap1 = candidateGap
                    break
            if centralGap1 == None:
                centralGap1 = gap1
                classFamily1.gapList.append(centralGap1)

            for candidateGap in classFamily2.gapList:
                if candidateGap.name == gap2.name:
                    centralGap2 = candidateGap
                    break
            if centralGap2 == None:
                centralGap2 = gap2
                classFamily2.gapList.append(centralGap2)

            #Now we know to which gap corresponds each gap found on the threads
            #It's time to rebuild the border references
            if 'left' in border:
                loggerMain.info('Left border found between {} and {}'.format(centralGap1.name, centralGap2.name))
                if centralGap2 not in centralGap1.leftMatch:
                    centralGap1.addBorderGapRef(centralGap2, 'left')
                if centralGap1 not in centralGap2.leftMatch:
                    centralGap2.addBorderGapRef(centralGap1, 'left')
            if 'right' in border:
                loggerMain.info('Right border found between {} and {}'.format(centralGap1.name, centralGap2.name))
                if centralGap2 not in centralGap1.rightMatch:
                    centralGap1.addBorderGapRef(centralGap2, 'right')
                if centralGap1 not in centralGap2.rightMatch:
                    centralGap2.addBorderGapRef(centralGap1, 'right')

    #To whoever cares: Code should never go into this loop, but no one knows for sure, so it will remain in order to
    # add some safety to the program. Its function is to remove all the gaps whose have no border gaps neither left
    # or right. Remove at your own risk
    for familyRoot in masterFamilyList:
        i = 0
        while i < len(familyRoot.gapList):
            gap  = familyRoot.gapList[i]
            if len(gap.leftMatch) + len(gap.rightMatch) == 0:
                familyRoot.gapList.remove(gap)
            else:
                i += 1

    saveHerringData('masterListTest.pk1')
    #loadedStuff = loadHerringData('masterListTest.pk1')
    #masterFamilyList = loadedStuff[0]
    #fileList = loadedStuff[1]

    checkForProblemType2(masterFamilyList, loggerMain)

    loggerMain.info('Final steps. Cleaning files created by the program')
    tempFastaFiles = ['*gap1left.temp.fasta', '*gap1right.temp.fasta', '*gap2left.temp.fasta', '*gap2right.temp.fasta',
                '*border1Reg.temp.fasta', '*border2Reg.temp.fasta', '*border1Adj.temp.fasta', '*border2Adj.temp.fasta']
    for file in os.listdir('.'):
        for pattern in tempFastaFiles:
            if fnmatch.fnmatch(file, pattern):
                os.remove(file)

    for file in os.listdir('./temp'):
        for pattern in tempFastaFiles:
            if fnmatch.fnmatch(file, pattern):
                os.remove(file)

    loggerMain.info('Storing all the info into two tipes of files')
    for family in masterFamilyList:
        family.writeGapInfo(family.findOrphans(type = ('left', 'right')), 'diagnostic-All' + '-' +
                            str(family.parents[0]) + '.txt')
        family.writeGapInfo([item for item in family.gapList if item in family.findOrphans(type=('both'))],
                            'diagnostic-Good' + '-' + str(family.parents[0]) + '.txt')

    #gapClusterList = checkForGapClusters(masterFamilyList)
    #dumpClusterData(gapClusterList)

    timeFinish = time.time()
    duration = timeFinish - timeStart
    if duration < 60:
        print('This program needed {} seconds'.format(duration))
    if duration > 60:
        print('This program needed {} hour(s) and {} minute(s) aprox'.format(duration//3600, round((duration%3600)/60)))
    loggerMain.info('This program needed {} hour(s), {} minute(s) and {} second(s)'.format(duration//3600,
                                                                (duration%3600)//60, round(((duration%3600)%60)/60)))
