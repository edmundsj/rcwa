import numpy as np

def complexNumberArrayFromString(stringRow):
    delimiter = 'i'
    rowOfStrings = stringRow.split(delimiter)
    rowOfStrings = [elem + "j" for elem in rowOfStrings]
    rowOfStrings.remove("j")
    rowOfStrings = np.array(rowOfStrings)
    rowOfComplexNumbers = rowOfStrings.astype(np.cdouble)

    return rowOfComplexNumbers

def numpyArrayFromFile(filename):
    """ Requires input file with all columns together on the same 18 rows """
    fileHandle = open(filename, 'r')
    delimiter = 'i'
    fileLines = fileHandle.readlines()
    data = None
    i = 0
    for line in fileLines:
        line = line.replace(" ", "")
        if line is not "":
            rowOfStrings = line.split(delimiter)
            rowOfStrings = [elem + "j" for elem in rowOfStrings]
            rowOfStrings.remove("\nj")
            rowOfStrings = np.array(rowOfStrings)
            rowOfComplexNumbers = rowOfStrings.astype(np.cdouble)
            if i is 0:
                data = rowOfComplexNumbers
            else:
                data = np.vstack((data, rowOfComplexNumbers))
            i += 1

    return data;

def numpyArrayFromSeparatedColumnsFile(filename):
    """ Requires an input file with columns 1 through 6 in the first 18 columns followed by a
    vertical spacer followed by columns 7 through 12 and so on """
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
    return data
