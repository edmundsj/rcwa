import numpy as np
from rcwa import Source


def assert_almost_equal(a, b, absoluteTolerance=1e-10, relativeTolerance=1e-9, errorMessage=""):
    if isinstance(a, Source):
        assert a == b
    np.testing.assert_allclose(a, b, atol=absoluteTolerance, rtol=relativeTolerance, err_msg=errorMessage);
    if isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
        assert(a.shape == b.shape)

    if np.isscalar(a):
        assert np.isscalar(a) and np.isscalar(b)
    else:
        np.testing.assert_equal(type(a), type(b))


def get_unequal_indices(a, b, absoluteTolerance=1e-10, relativeTolerance=1e-9):
    truthArray = np.greater(np.abs(a - b),
            (absoluteTolerance + relativeTolerance * np.abs(a - b)*np.ones(a.shape)))
    indexArray = np.argwhere(truthArray)
    return indexArray

