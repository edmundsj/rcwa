from rcwa import RectangularGrating, Crystal, TriangularGrating, Grating
import pytest
from numpy.testing import assert_allclose, assert_equal
import numpy as np

@pytest.fixture
def triangle():
    grating = TriangularGrating(thickness=0.2, er=2, er_void=0.5, ur=1, ur_void=1, Nx=3, Nz=3)
    return grating

def test_groove_width_nonsensical():
    with pytest.raises(ValueError):
        grating = RectangularGrating(groove_width=2, period=1)

def test_groove_data_er():
    grating = RectangularGrating()
    er_actual, ur_actual = grating._er_data(er=2, er_void=1, Nx=10, groove_fraction=0.6)
    er_desired= [1, 1, 1, 1, 1, 1, 2, 2, 2, 2]
    ur_desired = np.ones(10)

    assert_allclose(er_actual, er_desired)
    assert_allclose(ur_actual, ur_desired)

def test_groove_data_n():
    n = 2
    n_void = np.sqrt(2)
    grating = RectangularGrating()
    er_actual, ur_actual = grating._er_data(n=n, n_void=n_void, Nx=10, groove_fraction=0.349)

    er_desired= [2, 2, 2, 4, 4, 4, 4, 4, 4, 4]
    ur_desired = np.ones(10)

    assert_allclose(er_actual, er_desired)
    assert_allclose(ur_actual, ur_desired)

def test_groove_data_single():
    new_grating = RectangularGrating()
    n1 = 2
    n2 = 1

    material_array_desired = [2, 2, 2, 2, 2, 1, 1, 1, 1, 1]
    material_array_actual = new_grating._er_data_single(val1=n1, val2=n2, Nx=10, switch_fraction=0.5)
    assert_allclose(material_array_actual, material_array_desired)

def test_triangle():
    new_grating  = TriangularGrating(period=1, Nx=3, Nz=3)

def test_set_eun():
    grating = Grating()
    n, n_void, er, er_void, ur, ur_void = None, 0, 1, 2, 3, 4
    grating._set_eun(n=None, n_void=n_void, er=er, er_void=er_void, ur=ur, ur_void=ur_void)
    desired_vals = [np.sqrt(er*ur), np.sqrt(er_void*ur_void), er, er_void, ur, ur_void]
    actual_vals = [grating._n, grating._n_void, grating._er, grating._er_void, grating._ur, grating._ur_void]
    assert_allclose(actual_vals, desired_vals)

def test_set_eun_override():
    grating = Grating()
    n, n_void, er, er_void, ur, ur_void = 100, 9, 1, 2, 3, 4
    grating._set_eun(n=n, n_void=n_void, er=er, er_void=er_void, ur=ur, ur_void=ur_void)
    desired_vals = [n, n_void, np.square(n), np.square(n_void), 1, 1]
    actual_vals = [grating._n, grating._n_void, grating._er, grating._er_void, grating._ur, grating._ur_void]
    assert_allclose(actual_vals, desired_vals)

@pytest.mark.parametrize('period,lv,desired_period,desired_lv', [(1, None, 1, [1, 0]),(2, [1, 0], 1, [1, 0])])
def test_set_lattice_vector(period, lv, desired_period, desired_lv):
    grating = Grating()
    grating.set_lv_period(period=period, lattice_vector=lv)

    assert_equal(grating.period, desired_period)
    assert_equal(grating.lattice_vector, desired_lv)

def test_triangular_er_data(triangle):
    er_slices, ur_slices = triangle._er_data()
    desired_er_slices = [[0.5, 0.5, 0.5], [0.5, 0.5, 2],[0.5, 2, 2]]
    desired_ur_slices = [[1, 1, 1],[1, 1, 1], [1, 1, 1]]
    assert_equal(er_slices, desired_er_slices)
    assert_equal(ur_slices, desired_ur_slices)

def test_triangular_slice(triangle):
    layers = triangle.slice()
    desired_er_slices = [[0.5, 0.5, 0.5], [0.5, 0.5, 2],[0.5, 2, 2]]
    desired_ur_slices = [[1, 1, 1],[1, 1, 1], [1, 1, 1]]
    er_slices = [layer.crystal.permittivityCellData for layer in layers]
    ur_slices = [layer.crystal.permeabilityCellData for layer in layers]
    assert_equal(er_slices, desired_er_slices)
    assert_equal(ur_slices, desired_ur_slices)
