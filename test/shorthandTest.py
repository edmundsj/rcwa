import numpy as np

def assertAlmostEqual(a, b, absoluteTolerance=1e-10, relativeTolerance=1e-9, errorMessage=""):
    np.testing.assert_allclose(a, b, atol=absoluteTolerance, rtol=relativeTolerance, err_msg=errorMessage);

def assertEqual(a, b, errorMessage=""):
    np.testing.assert_equal(a, b, err_msg=errorMessage)

