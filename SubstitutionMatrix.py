"""
Creates a class that constructs a 2d array containing the substitution matrix
and has a method to return score by mapping
Avni Fatehpuria and Guillermo Barreto Corona, Spring 2017
"""

class SubstitutionMatrix(object):
    def __init__(self, inputFile):
        #constructs from a user input file name, no error checking since the method\
        #is called within a try/catch block

        matrix = open(inputFile, "r")
        valueString = matrix.readline()
        values = valueString.split()
        self.array = [[0 for x in range(len(values))] for y in range(len(values))]
        self.order = {}
        for i in range(len(values)):
        #the for i in len(values) is kinda a misnomer? we'll be using this to go down the
        #matrix and read in new lines, its just that that also has the same length
            scores = matrix.readline().split()
            self.order[scores[0]] = i
            for j in range(len(values)):
                self.array[i][j] = scores[j+1]
                #print self.array[i][j]


    def getScore(self, input1, input2):
        #first we look in order dict to know where in our array the two things we need are
        x = self.order[input1]
        y = self.order[input2]
        #then we index into that spot in the array and return what it stores
        return self.array[x][y]


def main():
    #just tests this code

    matrix = SubstitutionMatrix("blosum-62")
    print matrix.order

if __name__ == '__main__' :
    main()
