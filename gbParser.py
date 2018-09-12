from Bio import SeqIO
from itertools import filterfalse
import argparse
import sys
from subprocess import call



fosmidList = []

class Fosmid():
    def __init__(self, name, length, seq):
        self.seq = seq
        self.name = name
        self.length = length
        self.features = []
        self.featureDict = {}

    def addFeature(self, feature):
        #Check the feature type and add it to the appropiate list
        found = False
        for key in self.featureDict:

            if key == feature.type:
                self.featureDict[key] += 1
                feature.id = key.lower() + '_' + str(self.featureDict[key])
                found = True
                break
        #If the relevant type is not yet in the list, just add it
        if not found:
            self.featureDict[feature.type] = 1
            feature.id = feature.type + '_0'
        #Finally, add the feature to feature list
        self.features.append(feature)

    def purgeGeneList(self):
        #Some features appear 2 times in GenBank feature list: once as themselves and another as a gene. Both entries have the same
        # locus_tag qualifiers, so use those to remove the gene entries (CDS entries have more info)
        #In case this does not work correctly: db_xref is another possible qualifier, check both features' position

        locusList = []
        geneList = []
        otherList = []

        for feature in self.features:
            try:
                if feature.type != 'gene' and feature.qualifiers['locus_tag'] != None:
                    locusList.append(feature)
                else:
                    geneList.append(feature)
            except(KeyError):
                otherList.append(feature)

        print(len(locusList), len(geneList), len(otherList))

#add each feature whose doesnt exist in locusList
        locusList = sorted(locusList, key= lambda x:x.qualifiers['locus_tag'])
        newGeneList = [feature for feature in self.features if self._checkDuplicates(feature,locusList) == False]
        self.features = locusList + newGeneList

    def _checkDuplicates(self, feature, locuslist):
        #new meta
        pos_start = 0
        pos_end = len(locuslist)-1
        found = False
        while pos_start<pos_end and not found:
            try:
                pos_act = pos_start + (pos_end - pos_start) // 2
                #i += 1
                #print('inicio: {}, medio: {}, final: {}, iteracion: {}'.format(pos_start, pos_act, pos_end, i))
                if locuslist[pos_act].qualifiers['locus_tag'] == feature.qualifiers['locus_tag']:
                    found = True
                elif locuslist[pos_act].qualifiers['locus_tag'] < feature.qualifiers['locus_tag']:
                    if pos_start == pos_act:
                        break
                    pos_start = pos_act
                elif locuslist[pos_act].qualifiers['locus_tag'] > feature.qualifiers['locus_tag']:
                    pos_end = pos_act
            except(KeyError):
                return False
        return found

    def returnFeatureTypes(self):
        resultDict = {}
        for feature in self.features:
            if feature.type not in resultDict.keys():
                resultDict[feature.type] = 1
            else:
                resultDict[feature.type] += 1
        return resultDict

    def removeSourceFeature(self):
        for feature in self.features:
            if feature.type == 'source':
                self.features.remove(feature)
                break



class Feature():

    def __init__(self, Fosmid,  gbFeatList):
        self.id = None
        self.fosmid = Fosmid
        self.type = gbFeatList.type
        self.sequence = None
        #Get position
        self.position = [gbFeatList.location.start.position, gbFeatList.location.end.position, gbFeatList.location.strand]
        if self.position[2] == -1:
            self.position[2] = '-'
        elif self.position[2] == 1:
            self.position[2] = '+'
        self.qualifiers = gbFeatList.qualifiers


    def getFeatureSequence(self, sequence):
        self.sequence = sequence




#Get Filename
gbFiles = []
inputFiles = ['M1627.gbff']


'''
#Parse genbank,
for file in inputFiles:
    print(file)
    inputFile = SeqIO.parse(file, 'genbank')
    for record in inputFile:
        gbFiles.append(record)

print(len(gbFiles), 'records from', len(inputFiles), 'file(s) in input\n')

for gbRecord in gbFiles:

    newFosmid = Fosmid(name=gbRecord.id, length=gbRecord.features[0].location.end, seq=gbRecord.seq)
    featureList = gbRecord.features
    for rawFeature in featureList:
        newFeature = Feature(newFosmid, rawFeature)
        newFeature.getFeatureSequence(str(gbRecord.seq[rawFeature.location.start.position:rawFeature.location.end.position]))
        newFosmid.addFeature(newFeature)

    print(len(newFosmid.features))
    print(str(newFosmid.returnFeatureTypes()))

    newFosmid.purgeGeneList()

    print(str(newFosmid.returnFeatureTypes()))
'''

def getRecords(filenames):
    gbFiles = []
    for file in filenames:
        print(file)

        if tryNewFastaFile(file) is False:
            inputFile = SeqIO.parse(file, 'genbank')
            print('opening as genbank')
        else:
            inputFile = SeqIO.parse(file, 'fasta')
            print('opening as fasta')

        for record in inputFile:
            print(file,record.name, record.id, len(record.seq))
            gbFiles.append((file,record.name, record.id, len(record.seq)))
    return gbFiles

#new meta
def getNewRecords(filenames):
    gbFiles = []
    for file in filenames:
        name = None
        accession = None
        size = 0
        file_type = None
        i = 0
        added = False

        with open(file) as f:
            for line in f:
                if file_type == None:
                    content = line[:1]
                    if content == '>':
                        file_type = 'fasta'
                        print('opening as fasta')
                    elif content == 'L':
                        file_type = 'genBank'
                        print('opening as genbank')
                    else:
                        print("file extension not expected")
                        break
                if file_type == 'fasta':
                    if i == 0:
                        content = line.split()[0]
                        name = content[1:]
                        accession = content[1:]
                    else:
                        if line[:1] == '>':
                            print(file, name, accession, size)
                            gbFiles.append((file, name, accession, size))
                            content = line.split()[0]
                            name = content[1:]
                            accession = content[1:]
                            size = 0
                            i = 0
                        else:
                            # -1 because of the \n statement
                            size += len(line) - 1

                if file_type == 'genBank':
                    if i == 0:
                        content = line.split()
                        name = content[1]
                        size = content[2]
                    if i == 3:
                        content = line.split()
                        accession = content[1]
                    if not added and name != None and accession != None and size != 0:
                        print(file, name, accession, size)
                        gbFiles.append((file, name, accession, size))
                        added = True
                    elif added:
                        if len(line) == 2:
                            content = line.strip()
                            if content[0] == '/' and content[1] == '/':
                                added = False
                                i = -1
                i += 1
        if file_type == 'fasta':
            print(file, name, accession, size)
            gbFiles.append((file, name, accession, size))
    return gbFiles





def parseFastaFiles(filenames, exceptionDict = None):
    print('no of filenames:', len(filenames))
    fastaRecordList = []
    for file in filenames:
        with open(file, 'r') as filehandle:
            print('processing file')
            inputFile = SeqIO.parse(filehandle, 'fasta')
            for record in inputFile:
                print('processing record')
                for query in exceptionDict[file]:
                    if record.name == query[0] and len(record.seq) == int(query[2]):
                        fastaRecordList.append(record)
                    else:
                        print('match not found!')
                        print(record.name, query[0])
                        print(len(record.seq), int(query[2]))
    fosmidList = []
    print(len(fastaRecordList))
    for fastaRecord in fastaRecordList:
        newFosmid = Fosmid(name=fastaRecord.name, length=len(fastaRecord.seq), seq=fastaRecord.seq)
        fosmidList.append(newFosmid)
    return fosmidList


def parseGbFiles(filenames, exceptionDict = None):

    gbRecordList = []
    for file in filenames:
        with open(file, 'r') as filehandle:
            inputFile = SeqIO.parse(filehandle, 'genbank')
            for record in inputFile:
                for query in exceptionDict[file]:
                    if record.name == query[0] and record.id == query[1] and len(record.seq) == int(query[2]):
                        gbRecordList.append(record)
                    else:
                        print('match not found!')
                        print(record.name, query[0])
                        print(record.id, query[1])
                        print(len(record.seq), int(query[2]))
    fosmidList = []
    for gbRecord in gbRecordList:
        newFosmid = Fosmid(name=gbRecord.name, length=len(gbRecord.seq), seq=gbRecord.seq)
        featureList = gbRecord.features
        for rawFeature in featureList:

            newFeature = Feature(newFosmid, rawFeature)
            if (newFeature.position[1] - newFeature.position[0]) == newFosmid.length:
                continue
            else:
                newFeature.getFeatureSequence(
                    str(gbRecord.seq[rawFeature.location.start.position:rawFeature.location.end.position]))
                newFosmid.addFeature(newFeature)

        newFosmid.purgeGeneList()
        newFosmid.removeSourceFeature()
        fosmidList.append(newFosmid)

    return fosmidList

def tryNewFastaFile(filename):
    try:
        with open(filename) as f:
            content = f.readline()[:1]
            if content == '>':
                res = True
            else:
                res = False
        return res
    except Exception:
        return False




