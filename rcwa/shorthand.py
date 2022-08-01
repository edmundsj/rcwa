import numpy as np
import scipy as sp
import scipy.linalg
import math

inv = np.linalg.inv;
matrixExponentiate = sp.linalg.expm
matrixSquareRoot = sp.linalg.sqrtm
sqrt = np.lib.scimath.sqrt; # Takes sqrt of complex numbers successfully
sq = np.square;
eig = sp.linalg.eig # Performs eigendecomposition of identity intuitively (vectors are unit vectors)
norm = np.linalg.norm;
sin = np.sin;
cos = np.cos;
pi = np.pi;
dot = np.dot;
cross = np.cross;
diag = np.diag
diagonal = np.diagonal
conj = np.conj
real = np.real
imag = np.imag
deg = pi / 180
prod = np.prod

def fftn(data):
    """ Return the shifted version so the zeroth-order harmonic is in the center with
    energy-conserving normalization """
    dataShape = data.shape;
    return np.fft.fftshift(np.fft.fftn(data)) / np.prod(dataShape);

def complexArray(arrayInListForm):
    """ Wrapper for numpy array declaration that forces arrays to be complex doubles """
    return np.array(arrayInListForm, dtype=np.cdouble);

def complexIdentity(matrixSize):
    """ Wrapper for numpy identity declaration that forces arrays to be complex doubles """
    if matrixSize == 1:
        return 1
    else:
        return np.identity(matrixSize, dtype=np.cdouble);

def complexZeros(matrixDimensionsTuple):
    """ Wrapper for numpy zeros declaration that forces arrays to be complex doubles """
    return np.zeros(matrixDimensionsTuple, dtype=np.cdouble);

def complexOnes(matrixDimensionsTuple):
    return np.ones(matrixDimensionsTuple, dtype=np.cdouble);

def reshapeLowDimensionalData(data):
    if hasattr(data, '__len__'):
        dataShape = data.shape;
    else:
        dataShape = [1]
        data = np.array([data])

    if(len(dataShape) == 1): # we have only x-data.
        Nx = dataShape[0];
        data = data.reshape(Nx, 1, 1);
    elif(len(dataShape) == 2): # We have x and y data
            Nx = dataShape[0];
            Ny = dataShape[1];
            data = data.reshape(Nx, Ny, 1);
    elif(len(dataShape) == 3): # We have x- y- and z-data (
        data = data;
    else:
        raise ValueError(f"""Input data has too many ({len(dataShape)}) dimensions.
        Only designed for up to 3 spatial dimensions""");

    return data;

def kroneckerDeltaVector(size):
    vector = complexZeros(size)
    zeroLocation = math.floor(size/2)
    vector[zeroLocation] = 1
    return vector


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
        if line != "":
            rowOfStrings = line.split(delimiter)
            rowOfStrings = [elem + "j" for elem in rowOfStrings]
            rowOfStrings.remove("\nj")
            rowOfStrings = np.array(rowOfStrings)
            rowOfComplexNumbers = rowOfStrings.astype(np.cdouble)
            if i == 0:
                data = rowOfComplexNumbers
            else:
                data = np.vstack((data, rowOfComplexNumbers))
            i += 1

    fileHandle.close()
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
        if line != "":
            rowOfComplexNumbers = complexNumberArrayFromString(line)

            if rowNumber == 0:
                data[columnNumber] = rowOfComplexNumbers
            else:
                data[columnNumber] = np.vstack((data[columnNumber], rowOfComplexNumbers))
            rowNumber += 1

        if line == "": # This indicates we should start a new set of columns and append it to the old one
            columnNumber += 1
            rowNumber = 0

    fileHandle.close()

    data = np.hstack((data[0], data[1], data[2]))
    return data
