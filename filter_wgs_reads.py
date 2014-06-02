import sys, re
from time import time

if len(sys.argv) < 2:
    print 'usage: python', sys.argv[0], '<FILENAME>'
    print 'Script to filter pyrosequenced fasta files according to paper specifications'
    exit(0)

fileName = sys.argv[1]
filteredFileName = fileName.rsplit('.')[0] + '_filtered.fna'
# TODO: check file existence

inFile = open(fileName, 'r')
outFile = open(filteredFileName, 'w+')

seq = ''
head = ''

nFiltered = 0
nReads = 0
seqDict = {}

t0 = time()

def passesFilter(read):
    if len(read) < 60: return False
    degenerate = read.count('N')
    if degenerate > 2: return False
    if re.search('NN', read) != None: return False
    return True

# filter short and degenerate reads
for line in inFile:
    if line.startswith('>'):
        nReads += 1
        if passesFilter(seq):
            identifier = seq[0:20]
            seqDict.setdefault(identifier, []).append((head, seq))
        else:
            print 'filtered bad sequence'
            nFiltered += 1
        head = line
        seq = ''
    else:
        seq += line

print time() - t0, 'seconds'
print 'scanning', len(seqDict), 'possible duplicates'

# remove duplicates

def areDuplicates(readA, readB):
    sharedLength = min(len(readA), len(readB))
    identical = 0.0
    for i in xrange(sharedLength):
        if readA[i] == readB[i]:
            identical += 1.0
    return 1.0 * identical / sharedLength > 0.97

n = 0
for key in seqDict.keys():
    print 'step', n, time() - t0, 'seconds'
    n += 1
    candidates = seqDict[key]
    duplicates = set()
    if len(candidates) > 1:
        for i in xrange(len(candidates)-1):
            for j in xrange(i+1, len(candidates)):
                if areDuplicates(candidates[i][1], candidates[j][1]):
                    duplicates.add(candidates[j])
    nFiltered += len(duplicates)
    for dupe in duplicates:
        candidates.remove(dupe)

# print filtered fasta-file
for group in seqDict:
    for header, sequence in group:
        outFile.write(header + '\n')
        outFile.write(sequence + '\n')

print 'filtered', nFiltered, 'of', nReads, '=>', 100.0*(nReads-nFiltered)/nReads, 'percent kept'
