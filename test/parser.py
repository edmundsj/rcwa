import numpy as np

def complexParser(filename):
    fileHandle = open(filename, 'r')
    delimiter = ','
    fileLines = fileHandle.readlines()
    data = None
    i = 0
    for line in fileLines:
        line = line.replace(" ", "")
        line = line[:-2]; # remove the last trailing comma
        rowOfStrings = np.array(line.split(delimiter))
        rowOfComplexNumbers = rowOfStrings.astype(np.cdouble)
        if i is 0:
            data = rowOfComplexNumbers
        else:
            data = np.vstack((data, rowOfComplexNumbers))
        i += 1

    return data
