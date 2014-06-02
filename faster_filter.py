import sys, re

if len(sys.argv) < 2:
    print 'usage: python', sys.argv[0], '<FILENAME>'
    print 'Script to filter pyrosequenced fasta files according to paper specifications'
    exit(0)

fileName = sys.argv[1]
expr = re.compile('NN')

def getSeqCount(fileName):
    file = open(fileName, 'r')
    return len(file.read().split('>'))

def printPercentage(perc):
    global N_DROP, N_SEQS
    VERT = '0'
    HORIZ = '0'
    p = int(perc * 20)
    q = 20 - p
    string1 = sys.argv[0] + ' ' + sys.argv[1] + '\t- ' + str(N_SEQS) +'\n\n'
    string2 = ''
    try:
        string2 = '[' + '#'*p + '-'*q + '] ' + str(perc*100)[:5] + '%'
    except:
        pass
    string3 = '\t(' + str(N_DROP) + ' dropped)               \n\n'
    print('\033['+VERT+';'+HORIZ+'f' + string1 + string2 + string3)

def passesFilter(read):
    global expr
    if len(read) < 60: return False
    if read.count('N') > 2: return False
    if re.search('NN', read) != None: return False
    return True

def isUnique(read, theReads):
    reads = set(theReads)
    L = len(read)
    for elem in reads:
        readB = elem[1]
        identical = 0.0
        sharedLength = min(L, len(readB))
        for i in xrange(sharedLength):
            if read[i] == readB[i]:
                identical += 1.0
        if (identical / sharedLength) > 0.97:
            return False
    return True

def filterSeqs(fileName):
    global N_DROP, N_SEQS
    file = open(fileName, 'r')
    seq=''
    head=file.readline()
    N_READ = 0
    seqDict = {}
    for line in file:
        if line.startswith('>'):
            if not passesFilter(seq):
                N_DROP += 1
            else:
                key = hash(seq[:20])
                if key in seqDict:
                    if isUnique(seq, seqDict[key]):
                        seqDict[key].append((head, seq))
                        N_READ += 1
                    else:
                        N_DROP += 1
                else:
                    seqDict[key] = [(head, seq)]

            N_READ += 1
            if (N_READ % 100) == 0: printPercentage(1.0 * N_READ / N_SEQS)
            head = line
            seq = ''
        else:
            seq += line
    printPercentage(1.0)
    return seqDict

N_SEQS = getSeqCount(fileName)
N_DROP = 0

for n in xrange(25): print '\n'

filteredReads = filterSeqs(fileName)

outFileName = fileName.rsplit('.')[0] + '_filtered.fna'
outFile = open(outFileName, 'w+')
N_REMAINING = 0
for key in filteredReads.keys():
    read = filteredReads[key]
    for header, seq in read:
        N_REMAINING += 1
        outFile.write(header)
        outFile.write(seq)
print N_REMAINING, 'of', N_SEQS, '=', str(1.0*N_REMAINING/N_SEQS)[:6], 'percent remained\n'
