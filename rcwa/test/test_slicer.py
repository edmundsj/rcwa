from rcwa import Slicer
import pytest
from numpy.testing import assert_equal
import numpy as np

def test_valid_funcdata():
    with pytest.raises(ValueError) as e:
        Slicer()

def test_generate_coors():
    def test_func(x, y, z):
        return 1
    slicer = Slicer(func=test_func, Nx=2, Ny=2, Nz=2, xmin=0, xmax=1, ymin=0, ymax=1, zmin=0, zmax=1)
    x_desired = np.array([
       [[0, 0], [0, 0]],
        [[1, 1], [1, 1]]
    ])
    y_desired = np.array([
       [[0, 0], [1, 1]],
        [[0, 0], [1, 1]]
    ])
    z_desired = np.array([
       [[0, 1], [0, 1]],
        [[0, 1], [0, 1]]
    ])

    x, y, z = slicer.coordinates()
    assert_equal(x, x_desired)
    assert_equal(y, y_desired)
    assert_equal(z, z_desired)

def test_generate_coors_unequal():
    def test_func(x, y, z):
        return 1
    slicer = Slicer(func=test_func, Nx=2, Ny=2, Nz=2, xmin=0, xmax=1, ymin=0.1, ymax=0.9, zmin=0.2, zmax=0.8)
    x_desired = np.array([
       [[0, 0], [0, 0]],
        [[1, 1], [1, 1]]
    ])
    y_desired = np.array([
       [[0.1, 0.1], [0.9, 0.9]],
        [[0.1, 0.1], [0.9, 0.9]]
    ])
    z_desired = np.array([
       [[0.2, 0.8], [0.2, 0.8]],
        [[0.2, 0.8], [0.2, 0.8]]
    ])

    x, y, z = slicer.coordinates()
    assert_equal(x, x_desired)
    assert_equal(y, y_desired)
    assert_equal(z, z_desired)

def test_slice_square():
    def test_func(x, y, z):
        return np.square(x) + np.square(y) + np.square(z)
    slicer = Slicer(func=test_func, Nx=2, Ny=2, Nz=2, xmin=0, xmax=1, ymin=0, ymax=1, zmin=0, zmax=1)
    desired_vals = np.array([
       [[0, 1], [1, 2]],
        [[1, 2], [2, 3]]
    ])
    actual_vals = slicer.slice()
    assert_equal(actual_vals, desired_vals)
