import pytest
import numpy as np
from rcwa import Source
from rcwa.harmonics import generateKxMatrix, generateKyMatrix
from rcwa.testing import assertAlmostEqual

@pytest.fixture
def source():
    s = Source(theta=np.pi/2, phi=np.pi/6)
    yield s

def test_generate_kx(source):
    kx_desired = np.sqrt(3) / 2
    kx_actual = generateKxMatrix(source, None, 1)
    assertAlmostEqual(kx_desired, kx_actual)

def test_generate_ky(source):
    ky_desired = 1 / 2
    ky_actual = generateKyMatrix(source, None, 1)
    assertAlmostEqual(ky_desired, ky_actual)
