import numpy as np
import scipy as sp
import scipy.linalg
import math
from collections import OrderedDict

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
    return np.identity(matrixSize, dtype=np.cdouble);

def complexZeros(matrixDimensionsTuple):
    """ Wrapper for numpy zeros declaration that forces arrays to be complex doubles """
    return np.zeros(matrixDimensionsTuple, dtype=np.cdouble);

def complexOnes(matrixDimensionsTuple):
    return np.ones(matrixDimensionsTuple, dtype=np.cdouble);

def reshapeLowDimensionalData(data):
    dataShape = data.shape;
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
