from rcwa import Grating, Crystal
import pytest
from numpy.testing import assert_allclose
import numpy as np

def test_groove_width_nonsensical():
    with pytest.raises(ValueError):
        grating = Grating(groove_width=2, period=1)

def test_groove_data_er():
    grating = Grating()
    er_actual, ur_actual = grating._groove_data(er=2, er_void=1, Nx=10, groove_fraction=0.6)
    er_desired= [1, 1, 1, 1, 1, 1, 2, 2, 2, 2]
    ur_desired = np.ones(10)

    assert_allclose(er_actual, er_desired)
    assert_allclose(ur_actual, ur_desired)

def test_groove_data_n():
    n = 2
    n_void = np.sqrt(2)
    grating = Grating()
    er_actual, ur_actual = grating._groove_data(n=n, n_void=n_void, Nx=10, groove_fraction=0.349)

    er_desired= [2, 2, 2, 4, 4, 4, 4, 4, 4, 4]
    ur_desired = np.ones(10)

    assert_allclose(er_actual, er_desired)
    assert_allclose(ur_actual, ur_desired)

def test_groove_data_single():
    new_grating = Grating()
    n1 = 2
    n2 = 1

    material_array_desired = [2, 2, 2, 2, 2, 1, 1, 1, 1, 1]
    material_array_actual = new_grating._groove_data_single(val1=n1, val2=n2, Nx=10, switch_fraction=0.5)
    assert_allclose(material_array_actual, material_array_desired)

def test_grating_unknown_shape():
    with pytest.raises(ValueError) as e:
        new_grating = Grating(shape='nonsense')
