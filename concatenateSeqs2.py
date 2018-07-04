from Bio import SeqIO, Seq, SeqRecord, Alphabet
from subprocess import call
import os, sys, argparse, logging, time

#Initalize log
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG, filename='concatenateSeqs2-{}.log'.format(time.time()))


#Argparser stuff
argParser = argparse.ArgumentParser(description='concatenateSeqs joins a collection of contigs together using a complete sequence as reference. Requires blast 2.5+ and the Biopython library.\nWARNING: The algorithm used for the assembly is not perfect - if the output looks wrong, it\'s probably the assembler\'s fault.')
argParser.add_argument('-i', '--input', nargs = 1, required = True, type=str,
                       help='fasta file containing the contigs to be merged')
argParser.add_argument('-r', '--reference', nargs= 1, required = True, type=str,
                       help='fasta file containing the reference genome')
argParser.add_argument('-o', '--output', nargs=1, required = True, type=str,
                       help='prefix for the output files')
argParser.add_argument('--joinThreshold', nargs='?', default=500, type=int,
                       help='maximum distance between contigs to be joined together (default: 500). a value of -1 will join all contigs in order')
argParser.add_argument('--cutThreshold', nargs='?', default=0.02, type=float,
                       help= 'if >0, the script will attempt to cut contigs in order to align them to the reference. values between 0-1 (Default: 0.2)')
argParser.add_argument('--outfmt', nargs='?', default='fasta', type=str,
                       help='output format, either \'genbank\' or \'fasta\' (Default: fasta)')


#Initialize variables
blastpath = '/opt/ncbi-blast-2.6.0+/bin/'
refFasta = argParser.parse_args(sys.argv[1:]).reference[0]
contigFasta = argParser.parse_args(sys.argv[1:]).input[0]
outputPrefix = argParser.parse_args(sys.argv[1:]).output[0]
joinThreshold = argParser.parse_args(sys.argv[1:]).joinThreshold
cutThreshold = argParser.parse_args(sys.argv[1:]).cutThreshold
outputType = argParser.parse_args(sys.argv[1:]).outfmt
logger.info('Parameters passed to script: refFasta: {}, configFasta: {}, outputPrefix: {}, joinThreshold: {}, cutThreshold: {}, outFmt: {}'.format(refFasta, contigFasta,  outputPrefix,
            joinThreshold, cutThreshold, outputType))

contigList = {}
removedContigList = []
refList = []


def assembleContigs(contigList, threshold, niter = 1):

    assembleList = contigList

    for i in range(0, niter):
        logging.debug('\tstarting Contigs: {}'.format( len(assembleList)))
        assembleSingletons = []
        logging.debug('\tsingletons List: {}'.format(len(assembleSingletons)))
        assembleCandidates = []
        logging.debug('\tassembleCandidates: {}'.format(len(assembleCandidates)))
        contigList = sorted(assembleList, key=lambda contig: contig.refStart)
        for i in range(0, len(contigList) - 1):
            fstContig = contigList[i]
            scdContig = contigList[i + 1]
            try:

                if scdContig.refStart - fstContig.refEnd < threshold or threshold < 0:
                    logging.debug('\t{} {} {} {} - pass'.format(fstContig.id, fstContig.refEnd, scdContig.id, scdContig.refStart))
                    assembleCandidates.append([fstContig, scdContig])
                elif len(assembleCandidates) == 0:
                    logging.debug('\t{} not in list: {} < {}'.format(fstContig.id, scdContig.refStart - fstContig.refEnd, threshold ))
                    assembleSingletons.append(fstContig)

                elif fstContig not in assembleCandidates[-1]:
                    logging.debug(str(assembleCandidates.index(assembleCandidates[-1])))
                    logging.debug('\t{} not in {} , {}'.format(fstContig.id, assembleCandidates[-1][0].id,assembleCandidates[-1][1].id))
                    logging.debug('\t{} < {}'.format(scdContig.refStart - fstContig.refEnd, threshold))
                    assembleSingletons.append(fstContig)

            except (IndexError):
                print('Error!')





        # Assemble candidate pairs that have a contig in common
        logging.debug('\tn of assembly candidate pairs: {}'.format(len(assembleCandidates)))
        logging.debug('\tn of singletons:{}'.format(len(assembleSingletons)))
        if len(assembleCandidates) < 2:
            return assembleList
        else:


            i = 0
            while i < len(assembleCandidates) - 1:
                #print('i value', i)
                #print('number of pairs left', len(assembleCandidates))
                #print('attempting to merge pairs', i, i +1)
                #print(assembleCandidates[i + 1])
                if assembleCandidates[i][-1] == assembleCandidates[i + 1][0]:
                    #print('\tmerge!')
                    newList = assembleCandidates[i]
                    newList.append(assembleCandidates[i + 1][-1])
                    assembleCandidates[i] = newList
                    #print('\tnew list is', len(newList), 'seqs long, adding it to index',i)
                    assembleCandidates.remove(assembleCandidates[i + 1])
                    #print('removing list in index', i +1)
                    i = i
                    continue
                else:
                    i += 1
        mergedCandidates = []
        for list in assembleCandidates:
            newSeq = mergeContigs(list)
            mergedCandidates.append(newSeq)
        logging.debug('length of list: {}'.format(len(mergedCandidates)))
        del assembleList
        assembleList = mergedCandidates + assembleSingletons
        logging.info('length of list + singletons: {}'.format(len(assembleList)))
        del assembleSingletons
    return assembleList

def mergeContigs(list):
    contigList = sorted(list, key = lambda seq: seq.refStart)
    newRefSeq = contigList[0].refSeq
    newId = ''
    for i, contig in enumerate(contigList):
        if i < len(contigList) -1:
            newId += contig.id + '/'
        else:
            newId += contig.id
    newSeq = Seq.Seq('', alphabet=Alphabet.generic_dna)
    for contig in contigList:
        contig.seq.alphabet = Alphabet.generic_dna
        newSeq += contig.seq
    newRefStart = contigList[0].refStart
    newRefEnd = contigList[-1].refEnd

    newRecord = SeqRecord.SeqRecord(seq = newSeq, id = newId)
    newSeq = Sequence(newRecord, refSeq = newRefSeq)
    newSeq.refStart = newRefStart
    newSeq.refEnd = newRefEnd
    return newSeq

def getStrandCode(match, nextMatch):
    matchNo = 0
    nextMatchNo = 0
    if match.strand == 0:
        matchNo = -1
    else:
        matchNo = 0
    if nextMatch.strand == 0:
        nextMatchNo = 1
    else:
        nextMatchNo = 3
    return (matchNo + nextMatchNo)



class Sequence():

    def __init__(self, record, refSeq = None):
        self.refSeq = refSeq
        self.id = record.id
        self.seq = record.seq
        self.seq.description = 'None'
        self.seq.alphabet = Alphabet.generic_dna
        self.matchList = {}
        self.refStart = None
        self.refEnd = None
        self.strand = None
        self.length = len(self.seq)
        self.removedReason = None

    def addBlastHit(self, blastHit, refSeq):
        if refSeq.id not in self.matchList.keys():
            self.matchList[refSeq.id] = [blastHit]
        else:
            self.matchList[refSeq.id].append(blastHit)

    def processBlasts(self, blastFile, blastRef):
        for line in blastFile:
            splitLine = line.split('\t')
            if splitLine[0] in contigList.keys():
                contigList[splitLine[0]].addBlastHit(BlastHit(self, splitLine[3], splitLine[2], splitLine[8], splitLine[9], splitLine[6], splitLine[7], splitLine[11]), blastRef)

    def chooseRefSeq(self):
        logger.info('starting chooseRefSeq for contig {}'.format(self.id))
        candidate = None
        #If there is only blast for one  ref sequence, then there is no need for processing:
        if len(self.matchList.keys()) < 1:
            logger.warning('{} got no blasts, so it will be discarded'.format(self.id))
            self.removedReason = 'No blastHits assigned'
            self.refSeq = 'None'
            removedContigList.append(self)


        elif len(self.matchList.keys()) == 1:
            newList = []
            for key in self.matchList.keys():
                logger.info('only one blast, assigning {} as refSeq'.format(key))
                newList.append(key)
            self.refSeq = newList[0]

        #If there are a lot more blasts for one ref AND the avg bitscore is better, pick that one
        else:
            logger.info('more than one blast, starting bitScore calc')
            keyList = []
            for ref in self.matchList.keys():
                keyList.append(ref)
            keyDiv = len(self.matchList[keyList[0]]) / len(self.matchList[keyList[1]])
            keySubs = len(self.matchList[keyList[0]]) - len(self.matchList[keyList[1]])
            keyTotal = len(self.matchList[keyList[0]]) + len(self.matchList[keyList[1]])
            logger.debug('\tkeyDiv: {}, keySubs: {}, keyTotal:{}'.format(keyDiv,keySubs,keyTotal))

            bitScoreAvgList = []
            for ref in keyList:
                totalBitScore = 0

                for blastHit in self.matchList[ref]:
                    totalBitScore  += blastHit.bitScore
                totalBitScorefinal = totalBitScore / len(self.matchList[ref])
                logger.info('\ttotalBitScore for {}: {} / {} = {}'.format(ref, totalBitScore, len(self.matchList[ref]), totalBitScorefinal))
                bitScoreAvgList.append(totalBitScorefinal)

            if ((keyDiv > 2) or (keySubs > round(keyTotal * 0.20) and keyTotal > 20)) and (1.2 * bitScoreAvgList[0]) > bitScoreAvgList[1]:
                logger.info('\t{} selected as refSeq (bitScore + nº of blastHits)'.format(keyList[0]))
                self.refSeq = keyList[0]
            elif ((1/keyDiv) > 2 or (keySubs < round(keyTotal * 0.20) and keyTotal > 20)) and (1.2 * bitScoreAvgList[1]) > bitScoreAvgList[0]:
                logger.info('\t{} selected as refSeq (bitScore + nº of blastHits)'.format(keyList[1]))
                self.refSeq = keyList[1]
            #If both criteria can't be met, trust the bitscore
            elif (5 * bitScoreAvgList[0]) > bitScoreAvgList[1]:
                logger.info('\t{} selected as refSeq - based only on bitScore'.format(keyList[0]))
                self.refSeq = keyList[0]
            elif (5 * bitScoreAvgList[1]) > bitScoreAvgList[0]:
                logger.info('\t{} selected as refSeq - based only on bitScore'.format(keyList[1]))
                self.refSeq = keyList[1]
            else:
                logger.warning('{} was not assigned any refSeq and it will be discarded'.format(self.id))
                self.refSeq = 'None'
                self.removedReason = 'Not able to classify'
                removedContigList.append(self)


        #Remove the not chosen refSeq from the matchlist
        logger.info('\tremoving all refSeqs from matchlist that are not{}'.format(self.refSeq))
        if self.refSeq is not 'None':
            removeList = []
            for key in self.matchList.keys():
                if key is not self.refSeq:
                    removeList.append(key)
            for key in removeList:
                self.matchList.pop(key)

        if self.refSeq is not 'None':
            for seq in refList:
                if seq.id == self.refSeq:
                    self.refSeq = seq

    def filterBlasts(self, similarityThreshold, lengthThreshold, changeSelfList = False):
        'Keep matches that are at least 1% of the length ref'
        removedBlasts = 0
        newMatchList = {}
        for key in self.matchList.keys():
            newMatchList[key] = []
        for key in self.matchList.keys():
            for match in self.matchList[key]:
                if match.similarity >= similarityThreshold and match.lengthOfMatch > lengthThreshold:
                    newMatchList[key].append(match)
                else:
                    removedBlasts += 1
        if changeSelfList == True:
            self.matchList = newMatchList
        return newMatchList

    def findMaxEnd(self, matchList):
        '''Iterate through the matchlist and find the BlastHit that ends the furthest'''
        target = matchList[0].endRefPos
        queryTarget = 0
        for match in matchList:
            if match.endRefPos > target:
                target = match.endRefPos
                queryTarget = match.endQueryPos
        return [target, queryTarget]

    def findMinStart(self, matchList):
        '''Iterate through the matchlist and find the BlastHit that starts the earliest'''
        target = matchList[0].startRefPos
        queryTarget = 0
        for match in matchList:
            if match.startRefPos < target:
                target = match.startRefPos
                queryTarget = match.startQueryPos
        return [target, queryTarget]

    def writeFasta(self):
        with open('refSeq_temp.fasta', 'w') as newFastaFile:
            newFastaFile.write('>' + self.id + '\n')
            newFastaFile.write(str(self.seq))

    def groupBlasts(self):
        for keyRef in self.matchList:
            sortedRefBlasts = sorted(self.matchList[refSeq.id], key=lambda match: match.startRefPos, reverse=True)

            for i in range(0, len(sortedRefBlasts) -1):
                fstBlast = sortedRefBlasts[i]
                scdBlast = sortedRefBlasts[i + 1]

    def getOperatingMatchSize(self, startingThreshold, refSeq):
        threshold = startingThreshold
        logging.debug('{}, seqLength: {}'.format(self.id,self.length))
        #Bitscore includes length, so there is no need in analysing both
        rawList = self.filterBlasts(lengthThreshold=threshold, similarityThreshold=75.0)[refSeq.id]
        finalCandidates = []
        sortedBitScore = sorted(rawList, key=lambda match: match.bitScore, reverse=True)
        nOfCandidates = 0

        # Limit the candidates to the best 10 matches if there are more than 10
        logging.debug('\tn of total blasts: {}'.format(len(rawList)))
        if len(rawList) > 20:
            nOfCandidates = 10
        elif len(rawList) > 0:
            nOfCandidates = round(len(rawList) / 2) + 1
        logging.debug('\tnº of potential candidates: {}'.format(nOfCandidates))

        for i in range(0, nOfCandidates):
            # Get the best candidates of each category (nº depends on the number of matches)
            finalCandidates.append(sortedBitScore[i])
        logging.debug('\tnº of actual candidates: {}'.format(len(finalCandidates)))
        if len(finalCandidates) > 0:
            minStart = self.findMinStart(finalCandidates)[0]
            maxEnd = self.findMaxEnd(finalCandidates)[0]
            return([minStart, maxEnd])
        else:
            return None

    def placeContigOnRef(self, refSeq, startingThreshold, nOfIter):
        threshold = startingThreshold
        for i in range(0, nOfIter):
            logging.debug('\titer nº {}'.format(i+1))
            operatingSize = self.getOperatingMatchSize(threshold, refSeq)
            if operatingSize is not None:
                minStart = operatingSize[0]
                maxEnd = operatingSize[1]
                refSeqMin = len(refSeq.seq) * cutThreshold
                refSeqMax = len(refSeq.seq) - refSeqMin

                if (maxEnd-minStart) < int(2 * self.length):
                    self.refStart = minStart
                    self.refEnd = maxEnd
                    logging.debug('\tmaxEnd - minStart = {} < {}'.format(maxEnd-minStart, int(2 * self.length)))
                    logging.info('{} has found a good match in {}'.format(self.id, refSeq.id))
                    return(['OK', threshold])

                elif refSeqMin > minStart and refSeqMax < maxEnd:
                    logger.warning('{} is considered for cutting. refSeqMin ({}) > minStart ({}) & refSeqMax ({}) < maxEnd ({})'.format(self.id,
                                    refSeqMin, minStart, refSeqMax, maxEnd))
                    return(['CUT', threshold])

                else:
                    #logger.warning('{} is not considered for cutting. refSeqMin ({}) < minStart ({}) & refSeqMax ({}) > maxEnd ({})'.format(self.id, refSeqMin, minStart, refSeqMax, maxEnd))
                    threshold *= 1.5
                    logger.debug('\t{} - {}, calculatedLen: {}, seqLength: {}'.format(minStart, maxEnd, maxEnd - minStart, self.length))
                    logger.debug('\tCalculated contig length is larger than the actual contig, raising threshold to {}'.format(threshold))
                    continue

            else:
                    threshold /= 2
                    logger.debug('\tNo matches found, reducing threshold to {}'.format(threshold))
                    continue
        logger.warning('{} has reached max number of iterations and will be discarded'.format(self.id))
        return(['NO', threshold])

    def cutSequence(self, refSeq, threshold):
        logger.info('cutting sequence {}'.format(self.id))
        matchList = self.matchList[refSeq.id]
        matchList.sort(key=lambda match: match.startQueryPos)
        #Let's find relevant gaps between BlastHits!
        for i in range(0, len(matchList)-1):
            match = matchList[i]
            nextMatch = matchList[i + 1]
            # A gap between 2 blastHits is relevant when the dtce is almost as big as the ref in the reference sequence
            queryRefLengthDiff = abs(
                abs(nextMatch.startRefPos - match.endRefPos) - abs(nextMatch.startQueryPos - match.endQueryPos))
            if 0.9 < queryRefLengthDiff / refSeq.length < 1.1:
                logging.debug('\tfound a gap between matches {} and {}: abs(nextMatch.startRefPos({}) - match.endRefPos({})) - abs(nextMatch.startQueryPos({}) - match.endQueryPos({}))'.format(i,
                               i+1, nextMatch.startRefPos, match.endRefPos, nextMatch.startQueryPos, match.endQueryPos))
                logging.info('\t0.9 < queryRefLengthDiff ({}) / refSeq.length ({}) = {} > 1.1'.format(queryRefLengthDiff,  refSeq.length, queryRefLengthDiff / refSeq.length))
                # We found a gap, now let's check if it is the gap we're looking for
                #Maybe checkExtremeMatch is not needed
                #if match.checkExtremeMatch(refSeq) + nextMatch.checkExtremeMatch(refSeq) == 1:
                    # Matches are sorted by queryPos, so match will be always closer to the origin
                    #cutDtce = len(refSeq.seq) - (nextMatch.endQueryPos)

                cutDtce = 0
                strandCode = getStrandCode(match, nextMatch)

                if match.strand == 0:
                    dtceToEnd = refSeq.length - match.endRefPos
                    cutDtce = round((match.endQueryPos) + dtceToEnd)
                    logging.info('\tsequence cut! cutDtce = {} - {} = {}'.format(match.endQueryPos, dtceToEnd, cutDtce))
                else:
                    dtceToEnd = match.endRefPos
                    cutDtce = round((match.startQueryPos) + dtceToEnd)
                    logging.info('\tsequence cut! cutDtce = {} - {} = {}'.format(match.startQueryPos, dtceToEnd, cutDtce))
                #logging.info('\tsequence cut! cutDtce= (nextMatch.startQueryPos ({}) + match.endQueryPos ({})) / 2 = {}'.format(nextMatch.startQueryPos, match.endQueryPos, cutDtce))
                logging.debug('match.startQueryPos: {}, match.endQueryPos: {}, nextMatch.startQueryPos: {}, nextMatch.endQueryPos: {}'.format(
                            match.startQueryPos, match.endQueryPos, nextMatch.startQueryPos, nextMatch.endQueryPos))
                logging.debug('match.startRefPos: {}, match.endRefPos: {}, nextMatch.startRefPos: {}, nextMatch.endRefPos: {}'.format(
                        match.startRefPos, match.endRefPos, nextMatch.startRefPos, nextMatch.endRefPos))

                logging.debug('match strand: {}, nextMatch strand: {}'.format(match.strand, nextMatch.strand))

                operatingSize = self.getOperatingMatchSize(threshold, refSeq)

                seq1 = Sequence(SeqRecord.SeqRecord(self.seq[:cutDtce], id = self.id + '_1'), refSeq)

                seq2 = Sequence(SeqRecord.SeqRecord(self.seq[cutDtce:], id = self.id + '_2'), refSeq)
                if match.startRefPos < nextMatch.startRefPos:
                    seq1.refStart = 0
                    seq1.refEnd = seq1.length
                    seq2.refStart = refSeq.length - seq2.length
                    seq2.refEnd = refSeq.length
                else:
                    seq2.refStart = 0
                    seq2.refEnd = seq1.length
                    seq1.refStart = refSeq.length - seq2.length
                    seq1.refEnd = refSeq.length
                logging.info('2 new sequences: {} ({} - {}) w/ length {} and {} ({} - {}) w/ length {}'.format(seq1.id, seq1.refStart, seq1.refEnd, seq1.length, seq2.id,  seq2.refStart, seq2.refEnd, seq2.length))
                return([seq1, seq2])
            else:
                #logger.debug('\tcheckExtremeSequence returns: {} {}, going to the next pair of matchHits'.format(match.checkExtremeMatch(refSeq), nextMatch.checkExtremeMatch(refSeq)))
                logger.debug('gap is not of the correct size ({}), going to the next pair of matchHits'.format(queryRefLengthDiff / refSeq.length))


    def getStrand(self):
        strandNo = [0,0]
        for match in self.matchList:
            if match.strand == 0:
                strandNo[0] += 1
            else:
                strandNo[1] += 1
        if strandNo[0] > strandNo[1]:
            return 0
        else:
            return 1


class BlastHit():
    def __init__(self, parent, matchLen, similarity,  startPos, endPos, qStartPos, qEndPos, bitScore):
        self.parentSeq = parent
        self.lengthOfMatch = int(matchLen)
        self.similarity = float(similarity)
        self.startRefPos = int(startPos)
        self.endRefPos = int(endPos)
        self.startQueryPos = int(qStartPos)
        self.endQueryPos = int(qEndPos)
        self.strand = 0
        self.bitScore = bitScore

        #check the strand of the fragment, fix the start/end positions if it is from the complementary strand
        if self.startRefPos > self.endRefPos:
            self.strand = 1
            self.startRefPos = int(endPos)
            self.endRefPos = int(startPos)
        '''
        elif self.startQueryPos > self.endQueryPos:
            self.strand = 1
            self.startQueryPos = int(qEndPos)
            self.endQueryPos = int(qStartPos)
        '''


        #Fix bitscores in scientific notation

        if 'e' in self.bitScore:
            bitScoreSplit = self.bitScore.split('e')
            if bitScoreSplit[1][0] == '+':
                self.bitScore = float(bitScoreSplit[0]) * pow(10,  int(bitScoreSplit[-1]))

            elif bitScoreSplit[1][0] == '-':
                self.bitScore = float(bitScoreSplit[0]) * pow(10,  (int(bitScoreSplit[-1])* -1))

        else:
            self.bitScore = float(bitScore)

    def checkExtremeMatch(self, refSeq):
        logger.debug('starting checkExtremeMatch')
        refSeqMin = int(len(refSeq.seq) * cutThreshold)
        refSeqMax = len(refSeq.seq) - refSeqMin
        logger.debug('\tmin: {} < {} ({})'.format(self.startRefPos, refSeqMin, refSeq.length))
        logger.debug('\tmax: {} > {} ({})'.format(self.endRefPos, refSeqMax, refSeq.length))
        if self.startRefPos < refSeqMin:
            return 0
        elif self.endRefPos > refSeqMax:
            return 1
        else:
            return -1


#Get refSeqs & contigs

logger.info('loading ref file')
with open(refFasta, 'r') as refFile:
    parsedRefs = SeqIO.parse(refFile, 'fasta')
    for record in parsedRefs:
        refList.append(Sequence(record))
        logger.debug('\tref sequence {} loaded, with length {}'.format(record.id, len(record.seq)))

logger.info('loading fasta file')
with open(contigFasta, 'r') as contigFile:
    parsedContigs = SeqIO.parse(contigFile, 'fasta')
    for record in parsedContigs:
        contigList[record.id] = (Sequence(record))

#Perform blasts
logging.info('doing blasts')
for refSeq in refList:
    refSeq.writeFasta()
    call([blastpath + 'makeblastdb', '-in', 'refSeq_temp.fasta', '-out', 'dbTest', '-dbtype', 'nucl'])
    call([blastpath + 'blastn', '-query', contigFasta, '-db', 'dbTest', '-out', 'blastResults.blastn', '-num_threads',
          '4', '-outfmt', '6'])
    with open('blastResults.blastn', 'r') as blastFile:
        for contig in contigList:
            contigList[contig].processBlasts(blastFile, refSeq)
    os.remove('blastResults.blastn')
    os.remove('refSeq_temp.fasta')

#Check nº1
logger.info('number of reference sequences: {}'.format(len(refList)))
logger.info('number of contig sequences: {}'.format(len(contigList.keys())))
logger.info('choosing refSeqs for each contig')
for contig in contigList:
    contigList[contig].chooseRefSeq()



logger.info('choosing refSeq results:')
#Check nº2 (Post classifying seqs)
print(len(contigList.keys()))
countContigDict = {}
for contig in contigList:
    try:
        refSeq = contigList[contig].refSeq.id
    except AttributeError:
        refSeq = 'None'
    if refSeq not in countContigDict.keys():
        countContigDict[refSeq] = 1
    else:
        countContigDict[refSeq] += 1
for key in countContigDict.keys():
    logger.info('{}: {}'.format(key, countContigDict[key]))

#Create a new, clean contigList:
cleanContigList = []
for contig in contigList:
    if contigList[contig] not in removedContigList:
        cleanContigList.append(contigList[contig])
'''
for ref in refList:
    count = 0
    for contig in cleanContigList:
        if contig.refSeq == ref:
            count += 1
    print(ref.id, count)
'''


goodContigs = {}

logger.info('starting operatingMatchSize')

for i, ref in enumerate(refList):
    logger.info('starting operatingMatchSize for contigs assigned to {}'.format(ref.id))
    count = 0
    totalOKs = 0
    for contig in cleanContigList:
        if contig.refSeq == ref:
            logger.info('\tprocessing contig {} assigned to {}'.format(contig.id, contig.refSeq.id))
            count += 1
            check = contig.placeContigOnRef(contig.refSeq, startingThreshold= 3000, nOfIter= 10)
            if check[0] == 'OK':
                totalOKs += 1
                if ref.id not in goodContigs.keys():
                    goodContigs[ref.id] = [contig]
                else:
                    goodContigs[ref.id].append(contig)
            elif check[0] == 'NO':
                contig.removedReason = 'No relevant blast matches'
                removedContigList.append(contig)
            elif check[0] == 'CUT':
                newSeqs = contig.cutSequence(ref, check[1])
                if newSeqs is not None and cutThreshold > 0:

                    for seq in newSeqs:
                        totalOKs += 1
                        if ref.id not in goodContigs.keys():
                            goodContigs[ref.id] = [seq]
                            logging.info('Seq {} with refSeq {} added - new Dict'.format(seq.id, seq.refSeq.id))
                        else:
                            goodContigs[ref.id].append(seq)
                            logging.info('Seq {} with refSeq {} added'.format(seq.id, seq.refSeq.id))
                else:
                    logger.warning('{} cannot be cut, so it will be discarded'.format(contig.id))
                    contig.removedReason = 'Not a good match - Tried to cut'
                    removedContigList.append(contig)

    logger.info('results for {}: {}/{} placed. Total Contigs: {}'.format(ref.id, totalOKs, count,len(cleanContigList)))
    totalLength = 0

    for ref in refList:
        if ref.id not in goodContigs.keys():
            goodContigs[ref.id] = []

    for contig in goodContigs[ref.id]:
        totalLength += contig.length


    logger.info('assemble completion: {}%'.format((totalLength / ref.length)* 100))


logger.info('')
logger.info('')


assembledContigs = {}
for ref in refList:
    logger.info('start assembling contigs for {}'.format(ref.id))
    assembledContigs[ref.id] = assembleContigs(goodContigs[ref.id], joinThreshold)

for ref in refList:
    if ref.id not in assembledContigs.keys():
        assembledContigs[ref.id] = []

for ref in assembledContigs:
    logger.info('{}: {}'.format(ref, len(assembledContigs[ref])))
    nuList = assembledContigs[ref]
    nuList = sorted(nuList, key = lambda contig: len(contig.seq), reverse=True)
    for contig in nuList:
        logger.info('\t {}: {}'.format(contig.id, contig.length))

logger.info('writing output Files')
for ref in assembledContigs:
    goodList = []
    for i, contig in enumerate(assembledContigs[ref]):
        goodList.append(SeqRecord.SeqRecord(contig.seq, id='contig_'+ str(i)))
    with open(outputPrefix + ref + '.goodContigs.' + outputType, 'w') as outputGood:
        SeqIO.write(goodList, outputGood, outputType)


badList = []
for contig in removedContigList:
    badList.append(SeqRecord.SeqRecord(contig.seq, contig.id))
with open(outputPrefix + '.badContigs.' + outputType, 'w') as outputBad:
    SeqIO.write(badList, outputBad, outputType)











