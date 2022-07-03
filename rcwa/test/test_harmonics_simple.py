import pytest
import numpy as np
from rcwa import Source, Crystal
from rcwa.harmonics import kx_matrix, ky_matrix, x_components, y_components
from rcwa.testing import assert_almost_equal

@pytest.fixture
def source():
    s = Source(theta=np.pi/2, phi=np.pi/6)
    yield s

@pytest.fixture
def crystal_1D():
    er = np.array([1, 1, 1, 3, 3, 3])
    ur = 1 + 0 * er
    lattice_vector = [2, 0]
    crystal = Crystal(lattice_vector, er=er, ur=ur)

    yield crystal

def test_generate_kx_0D(source):
    kx_desired = np.sqrt(3) / 2
    kx_actual = kx_matrix(source, None, 1)
    assert_almost_equal(kx_desired, kx_actual)

def test_generate_ky_0D(source):
    ky_desired = 1 / 2
    ky_actual = ky_matrix(source, None, 1)
    assert_almost_equal(ky_desired, ky_actual)

def test_generate_kx_1D(source, crystal_1D):
    kx_incident = x_components(source.k_incident)
    reciprocal_vector_x = x_components(crystal_1D.reciprocal_lattice_vectors[0])
    Kx_desired = np.diag(kx_incident - reciprocal_vector_x * np.array([-1, 0, 1]))
    Kx_actual = kx_matrix(source, crystal_1D, 3)

    assert_almost_equal(Kx_desired, Kx_actual)
