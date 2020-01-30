import numpy as np

def complexNumberArrayFromString(stringRow):
    delimiter = 'i'
    rowOfStrings = stringRow.split(delimiter)
    rowOfStrings = [elem + "j" for elem in rowOfStrings]
    rowOfStrings.remove("j")
    rowOfStrings = np.array(rowOfStrings)
    rowOfComplexNumbers = rowOfStrings.astype(np.cdouble)

    return rowOfComplexNumbers

# This is not working. For some reason data is being concatenated into a single array after 
# I assign the zeroth element to be a numpy array. I don't know how to fix this.
# This is definitely a problem with python or some weird subtle bug. Still need to fix.
def numpyArrayFromSeparatedColumnsFile(filename):
    fileHandle = open(filename, 'r')
    fileLines = fileHandle.readlines()
    data = [None, None, None]
    rowNumber = 0
    columnNumber = 0

    for line in fileLines:
        line = line.replace(" ", "")
        line = line.replace("\n", "")
        if line is not "":
            rowOfComplexNumbers = complexNumberArrayFromString(line)

            if rowNumber is 0:
                data[columnNumber] = rowOfComplexNumbers
            else:
                data[columnNumber] = np.vstack((data[columnNumber], rowOfComplexNumbers))
            rowNumber += 1

        if line is "": # This indicates we should start a new set of columns and append it to the old one
            columnNumber += 1
            rowNumber = 0

    data = np.hstack((data[0], data[1], data[2]))
    return np.array(data)

data = complexParser("test/matrixData/S21Layer1.txt")
print(data.shape)
