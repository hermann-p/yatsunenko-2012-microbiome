import os

biomNames = os.listdir('.')

scriptFile = open('combine.sh', 'w+')
scriptFile.write('#!/bin/bash\n')

firstRun = True
while len(biomNames) > 0:
    thisRun = []
    string = 'merge_otu_tables.py -o combined.biom -i '
    if not firstRun:          # concat combined.biom only when not empty
        string += 'combined.biom,'
    else:
        firstRun = False
    for n in range(50):
        name = None
        try:
            name = biomNames.pop()
            if name == 'combined.biom' or not name.endswith('.biom'):
                name = None   # filter undesired input files
        except:               # catch pop()ing from empty list
            pass
        if name != None:
            thisRun.append(name)
    string += ','.join(thisRun) + '\n' # append input files as comma-seperated list
    scriptFile.write(string)
