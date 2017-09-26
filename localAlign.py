"""
Avni Fatehpuria and Guillermo Barreto Corona
The purpose of this program is take in 2 sequences, a gap penalty, and a
filename and output an alignment or alignments for those 2 sequences based on
the given gap penalty. the file also creates a file with the raw sequences
(gaps included) of the alignements.
"""
from SubstitutionMatrix import *
from sys import *

"""
 usage - prints out information for how program should be used on the command
    line.  Call this method if the user gives improper arguments
"""
def main():
    inputs = getInput(argv)
    if (inputs == "none"):
        return
    seqList = printAlignment(inputs)
    makeFile(seqList,inputs[4])

def usage():
    print >> stderr, "Usage: python localAlign.py seq1 seq2 matrixFile g output"
    print >> stderr, "  seq1, seq2 - fasta files containing two sequences to align (1 per file)"
    print >> stderr, "  matrixFile - file containing substitution matrix"
    print >> stderr, "  g - integer specifying penalty for a gap"
    print >> stderr, "  output - name of output file for raw alignment"

#TODO: define other helper methods and a main method below

def getInput(inpt):
    #"parses the user input into usable things, returns a list of them"

    if (len(inpt)!=6):
         print "wrong number of arguments given"
         usage()
         return "none"

    try:
        sequence1 = ""
        fasta1 = open(inpt[1], "r")
        fasta1.readline()
        for line in fasta1:
            line = line.strip()
            sequence1 += str(line)
    except:
        print "failed to read in the first fasta file"
        usage()
        return "none"

    try:
        sequence2 = ""
        fasta2 = open(inpt[2], "r")
        fasta2.readline()
        for line in fasta2:
            line = line.strip()
            sequence2 += str(line)
    except:
        print "failed to read in the second fasta file"
        usage()
        return "none"

    try:
        subMatrix = SubstitutionMatrix(inpt[3])
    except:
        print "failed to read in the substitution Matrix"
        usage()
        return "none"

    try:
        gap = float(inpt[4])
    except:
        print "failed to read in the gap penalty"
        usage()
        return "none"

    try:
        output = inpt[5]
    except:
        print "no output file"
        usage()
        return "none"

    return ([sequence1, sequence2, subMatrix, gap, output])

#scoredict = {}
prev = {}

def findScore(i, j, sub, gap, scoreMatrix):
    """finds the best possible score for a single spot in the matrix
    uses highroad protocol to break ties"""
    diagonal = scoreMatrix[i-1][j - 1] + int(sub)
    left = scoreMatrix[i][j-1] + gap
    vertical = scoreMatrix[i-1][j] + gap
    if (max(0, diagonal, vertical, left) == 0):
        prev[(i, j)] = "none"
        return 0
    if (max(diagonal, vertical, left) == vertical):
        prev[(i, j)] = (i-1, j)
        return vertical
    if (max(diagonal, left) == diagonal):
        prev[(i, j)] = (i-1,j-1)
        return diagonal
    prev[(i, j)] = (i, j-1)
    return left

def scoring_matrix(seq1, seq2, subMatrix, gap):
    """ loops through a matrix to find the best alignment score of the sequences"""
    scoreMatrix = [[0 for x in range(len(seq2)+ 1)] for y in range(len(seq1)+1)]
    maxScore = 0
    maxLocation = [[0, 0]]
    for i in range(len(seq1)+ 1):
        scoreMatrix[i][0] = 0
        prev[(i, 0)] = "none"
        #scoredict[i, 0] = 0
    for j in range(len(seq2)+1):
        scoreMatrix[0][j] = 0
        prev[(0,j)] = "none"
        #scoredict[0,j] = 0
    for i in range(1, len(seq1)+ 1):
        for j in range(1, len(seq2)+1):
            scoreMatrix[i][j] = findScore(i, j, subMatrix.getScore(seq1[i-1], seq2[j-1]), gap, scoreMatrix)
            if (scoreMatrix[i][j] > maxScore):
                maxScore = scoreMatrix[i][j]
                maxLocation = [[i, j]]
            elif (scoreMatrix[i][j] == maxScore):
                maxLocation.append([i, j])

    return scoreMatrix, maxLocation

def traceback(scores, maxLocation, seq1, seq2):
    """
    This function traces through the matrix and returns the optimal
    alignments.
    """
    i = maxLocation[0]
    j = maxLocation[1]
    num1 = i
    num2 = j
    seq1result = ""
    seq2result = ""
    previous = prev[(i, j)]
    while (previous!="none"):
        if (previous[0] == i-1 and previous[1] == j-1):
            seq1result = seq1[i-1] + seq1result
            seq2result = seq2[j-1] + seq2result
            j = j - 1
            i = i - 1
            num1 = num1 -1
            num2 = num2 -1
        elif (previous[0] == i-1 and previous[1] == j):
            seq2result = "-" + seq2result
            seq1result = seq1[i-1] + seq1result
            i = i - 1
            num1 = num1 -1
        else:
            seq1result = "-" + seq1result
            seq2result = seq2[j-1] + seq2result
            j = j-1
            num2 = num2 -1
        previous = prev[(i,j)]

    return seq1result, seq2result, num1+1, num2+1

def makeFile(seqList, filename):
    """
    This function takes in a sequenceList and a filename. The sequenceList is a
    list of sequence alignments we want to write on the file. We create a file
    with the given filename.
    """
    f = open(str(filename), "w+")
    length = len(seqList)
    for x in range(0, length, 2):
        f.write(seqList[x]+"\n")
        f.write(seqList[x+1]+"\n")
        f.write("\n")


def printAlignment(inputs):
    """
    This function takes in inputs from argv and returns a list of sequences
    alignments. The purpose of this function, however, is to print the
    alignments in a readable, easy to look at manner.
    """
    seq1 = inputs[0]
    seq2 = inputs[1]
    matrix = inputs[2]
    gap = inputs[3]

    #Keeps a list of the ouput sequences
    seqList = []

    scores, maxLocation = scoring_matrix(seq1, seq2, matrix, gap)

    print "______________________________________________________________________________"

    for item in maxLocation:
        seqres1, seqres2, j, i = traceback(scores, item, seq1, seq2)
        seqList.append(seqres1)
        seqList.append(seqres2)
        string = ""

        #length of alignment
        leng = len(seqres2)

        #Accumalator that keeps track of num identities and gaps
        identity = 0
        gap = 0

        #Creates string with "|" and "." that correspond to the alignment
        for x in range(len(seqres1)):
            if seqres1[x] == seqres2[x]:
                string += "|"
                identity +=1
            elif seqres1[x] == "-" or seqres2[x] == "-":
                string += " "
                gap += 1
            else:
                string += "."
        #Stylistic choices ~
        print ""
        print "Alignment Score: " + str(int(scores[item[0]][item[1]]))
        print "______________________________________________________________________________"
        print ""
        print "Alignment:"
        print ""

        #blank1 and blank2 keep track of num gaps per 50 indices
        blank1 = 0
        blank2 = 0

        #start and end keep track of indices within the sequence
        start = 0
        end = 49

        while leng > 0:
            #We want to make sure that 50 indices in each line
            if leng >= 50:
                #reset blanks
                blank1 = 0
                blank2 = 0
                #keeps track of num gaps
                for x in range(1, 50):
                    if seqres1[start+x] == "-":
                        blank1 += 1
                    if seqres2[start+x] == "-":
                        blank2 += 1
                #Following code makes sure alignments print 'pretty'
                space1 = ""
                space2 = ""
                space3 = ""
                if len(str(j)) > len(str(i)):
                    for x in range(len(str(j))):
                        space2 += " "
                    for x in range(len(str(j))-len(str(i))):
                        space3 += " "
                elif len(str(i)) > len(str(j)):
                    for x in range(len(str(i))):
                        space2 += " "
                    for x in range(len(str(i))-len(str(j))):
                        space1 += " "
                else:
                    for x in range(len(str(i))):
                        space2 += " "
                #Gets end index for each line
                endl1 = j + (50-blank1)
                endl2 = i + (50-blank2)

                print str(j) + space1 + "  " + seqres1[start:end] + "  " + str(endl1)
                print "  " + space2 + string[start:end]
                print str(i) + space3 + "  " + seqres2[start:end] + "  " + str(endl2)
                print ""
                print ""

                #increases indices, redefines start, end, decreases leng
                j += (50-blank1)
                i += (50-blank2)
                start += 50
                end += 50
                leng = leng - 50

            #else statement account when less than 50 amino acids left
            #code is super similar to as above
            else:
                blank1 = 0
                blank2 = 0
                for x in range(1, len(seqres1[start:-1])):
                    if seqres1[start+x] == "-":
                        blank1 += 1
                for x in range(1, len(seqres2[start:-1])):
                    if seqres2[start+x] == "-":
                        blank2 += 1
                space1 = ""
                space2 = ""
                space3 = ""
                if len(str(j)) > len(str(i)):
                    for x in range(len(str(j))):
                        space2 += " "
                    for x in range(len(str(j))-len(str(i))):
                        space3 += " "
                elif len(str(i)) > len(str(j)):
                    for x in range(len(str(i))):
                        space2 += " "
                    for x in range(len(str(i))-len(str(j))):
                        space1 += " "
                else:
                    for x in range(len(str(i))):
                        space2 += " "
                len1 = len(seqres1[start:-1])
                len2 = len(seqres2[start:-1])
                endl1 = j + len1 - blank1
                endl2 = i + len2 - blank2
                print str(j) + space1 + "  " + seqres1[start:-1] + "  " + str(endl1)
                print space2 + "  " + string[start:-1]
                print str(i) + "  " + space3 + seqres2[start:-1] + "  " + str(endl2)
                leng = 0
        print ""
        print ""
        #prints identities and gaps
        print "Identity " + str(identity) + "/" + str(len(seqres1)) + " " + "Gaps " + str(gap) + "/" + str(len(seqres1))
        print "______________________________________________________________________________"
    return seqList



main()
