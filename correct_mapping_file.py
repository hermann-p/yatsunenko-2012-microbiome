import os, re
from fnmatch import fnmatch

INFILENAME='/home/pah11816/modified_decimaldot.csv'
OUTFILENAME = '/home/pah11816/Yatsunenko/Data/mapping.csv'
PATH_TO_READS='/home/pah11816/Yatsunenko/Data/16S'

def getTable(inFileName, outFile):
    theTable = {}
    file = open(inFileName, 'r')
    for line in file:
        if line.startswith('#'):
            outFile.write(line)
            continue
        try:
            stuff = line.split('\t')
            key = stuff[0]
            value = '\t'.join(stuff[1:])
            theTable[key] = value
        except IndexError:
            pass
    return theTable

def modifyIDs(table, theOutFile):
    notFound = 0
    outFile = theOutFile
    for root, dirs, files in os.walk(PATH_TO_READS):
        for basename in files:
            if (fnmatch(basename, '*.fna')):
                fileName = os.path.join(root, basename)
                tmpFile = open(fileName, 'r')
                # format: ">idString.idNum_readNum", so split it as follows
                longID = tmpFile.readline()[1:].split('_')[0]
                elements = longID.split('.')
                key = '.'.join(elements[:-1])
                idNum = elements[-1]
                mapEntry = ''
                try:
                    mapEntry = key + '.' + idNum + '\t' + table[key]
                    outFile.write(mapEntry)
                except KeyError: # table is not 100% consistent with read names
                    key = '.'.join(elements[:-2])
                    idNum = '.'.join(elements[-2])
                    try:
                        mapEntry = key + '.' + idNum + '\t' + table[key]
                        outFile.write(mapEntry)
                    except KeyError: # 7 samples were excluded from analysis
                        print key, '(' + longID + ') not found, giving up'
                        notFound += 1
    print notFound, 'keys not found'

outFile = open(OUTFILENAME, 'w+')
table = getTable(INFILENAME, outFile)
modifyIDs(table, outFile)