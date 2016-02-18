# Written by Alexander M. West
# Version: Dec. 4, 2015

#The consensus maker uses a simple "majority/plurality rules" algorithm to make a consensus at each base position.

def consensusMaker (file):
    nucKeyDict = {0:'T', 1:'C', 2:'G', 3:'A', 4:'-'}
    nucIdentityList = {}
    totalPositions = 0
    consensusRead = ''
    currentSeq = ""
    #agreementThreshold = float (5)/10  #require 50% agreement to call a bp - not currently in use

    f = open(file)
    for line in f:
        #print currentSeq
        if '>' in line and currentSeq != '':
            if totalPositions == 0:
                totalPositions = len(currentSeq) - 1
                for i in range(totalPositions):
                    nucIdentityList[i] = [0, 0, 0, 0, 0, 0] # In the order of T, C, G, A, -, Total
            for character in range(len(currentSeq)):
    #for i in xrange(readLength) : # Count the types of nucleotides at a position in a read. i is the nucleotide index within a read in groupedReadsList
    #for j in xrange(len(groupedReadsList)): # Do this for every read that comprises a SMI group. j is the read index within groupedReadsList
                try:
                    if currentSeq[character] == 't' :
                        nucIdentityList[character][0] += 1
                    elif currentSeq[character] == 'c':
                        nucIdentityList[character][1] += 1
                    elif currentSeq[character] == 'g':
                        nucIdentityList[character][2] += 1
                    elif currentSeq[character] == 'a':
                        nucIdentityList[character][3] += 1
                    elif currentSeq[character] == '-':
                        nucIdentityList[character][4] += 1
                    nucIdentityList[character][5] += 1
                except:
                    break
            currentSeq = ""
        else:
            currentSeq += line.rsplit()[0]
    #tally the last sequence, avoid fencepost issue
    for character in range(len(currentSeq)):
        try:
            if currentSeq[character] == 't' :
                nucIdentityList[character][0] += 1
            elif currentSeq[character] == 'c':
                nucIdentityList[character][1] += 1
            elif currentSeq[character] == 'g':
                nucIdentityList[character][2] += 1
            elif currentSeq[character] == 'a':
                nucIdentityList[character][3] += 1
            elif currentSeq[character] == '-':
                nucIdentityList[character][4] += 1
            nucIdentityList[character][5] += 1
        except:
            break
    for i in range(totalPositions):
        #print nucIdentityList[i]
        try:
            #Order of testing is -, G, C, A, T to give bias to earlier items in list.
            pluralityRead = 4
            pluralityQuantity = int(nucIdentityList[i][4])
            for j in [2, 1, 3, 0] :
                if nucIdentityList[i][j] > pluralityQuantity:
                    pluralityRead = j
                    pluralityQuantity = nucIdentityList[i][j]
            consensusRead += nucKeyDict[pluralityRead]
                #print float(nucIdentityList[i][j])/float(nucIdentityList[i][5])
                #if float(nucIdentityList[i][j])/float(nucIdentityList[i][5]) >= agreementThreshold :
                #consensusRead += nucKeyDict[j]
                #break
                #elif j==0:
                #if float(nucIdentityList[i][2] + nucIdentityList[i][3])/float(nucIdentityList[i][5]) >= agreementThreshold:
                #    consensusRead += "R"
                #        break
                #            elif float(nucIdentityList[i][1] + nucIdentityList[i][0])/float(nucIdentityList[i][5]) >= agreementThreshold:
                #        consensusRead += "Y"
                #        break
                #    elif float(nucIdentityList[i][1] + nucIdentityList[i][2])/float(nucIdentityList[i][5]) >= agreementThreshold:
                #        consensusRead += "S"
                #        break
                #    elif float(nucIdentityList[i][0] + nucIdentityList[i][3])/float(nucIdentityList[i][5]) >= agreementThreshold:
                #        consensusRead += "W"
                #        break
                #    elif float(nucIdentityList[i][0] + nucIdentityList[i][2])/float(nucIdentityList[i][5]) >= agreementThreshold:
                #        consensusRead += "K"
                #        break
                #    elif float(nucIdentityList[i][3] + nucIdentityList[i][1])/float(nucIdentityList[i][5]) >= agreementThreshold:
                #        consensusRead += "M"
                #        break
                #    else:
                #        consensusRead += "N"
                #        break
        except:
            consensusRead += "N"
    return consensusRead