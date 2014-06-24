def checkSampleIDmapping(fileName):
    print "checking", fileName, "-->",
    csv = open(fileName, "r")
    count = 0
    errors = 0
    csv.readline() # skip header

    for line in csv:
        line_split = line.split('\t')
        try:
            longID = line_split[0]
            shortID = line_split[1]
            if not longID.startswith(shortID):
                print "error in line", count, shortID, longID
                errors += 1
        except:
            pass
        count += 1

    print "read", count, "samples,", errors, "errors"

for fileName in ["16s_mapping_decimaldot.csv", "wgs_mapping_decimaldot.csv"]:
    checkSampleIDmapping(fileName)
