import numpy as np
import sys
sys.path.append('core')
from source import Source

def assertAlmostEqual(a, b, absoluteTolerance=1e-10, relativeTolerance=1e-9, errorMessage=""):
    if isinstance(a, Source):
        a = a.getRepresentationVector()
        b = b.getRepresentationVector()
    np.testing.assert_allclose(a, b, atol=absoluteTolerance, rtol=relativeTolerance, err_msg=errorMessage);
    if isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
        assert(a.shape == b.shape)

def getUnequalIndices(a, b, absoluteTolerance=1e-10, relativeTolerance=1e-9):
    truthArray = np.greater(np.abs(a - b),
            (absoluteTolerance + relativeTolerance * np.abs(a - b)*np.ones(a.shape)))
    indexArray = np.argwhere(truthArray)
    return indexArray

def assertEqual(a, b, errorMessage=""):
    np.testing.assert_equal(a, b, err_msg=errorMessage)

def assertArrayEqual(a, b, errorMessage=""):
    np.testing.assert_array_equal(a, b, err_msg=errorMessage)

def assertStringEqual(a, b):
    np.testing.assert_array_equal(a, b)
