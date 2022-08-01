from rcwa import Source, Crystal
from rcwa.shorthand import *
from rcwa.testing import *
from numpy.testing import assert_equal, assert_allclose
import pytest

@pytest.fixture
def crystal_1D():
    t1 = complexArray([1, 0])
    square_crystal = Crystal(t1)
    return square_crystal

@pytest.fixture()
def square_3d():
    t1  =  complexArray([1, 0, 0])
    t2 = complexArray([0,1,0])
    t3 = complexArray([0, 0, 1])
    oblique_crystal = Crystal(t1, t2, t3)
    yield oblique_crystal


@pytest.fixture()
def rectangular_3d():
    t1 = complexArray([1, 0, 0])
    t2 = complexArray([0,2,0])
    t3 = complexArray([0, 0, 1])
    oblique_crystal = Crystal(t1, t2, t3)
    yield oblique_crystal


@pytest.fixture()
def oblique_3d():
    t1  =  complexArray([1, 0, 0])
    t2 = complexArray([1,2,1])
    t3 = complexArray([2, 1, 1])
    oblique_crystal = Crystal(t1, t2, t3)
    yield oblique_crystal

@pytest.fixture
def oblique_2d():
    t1 = complexArray([0, 1])
    t2 = complexArray([1, 1])
    oblique_crystal = Crystal(t1, t2)
    yield oblique_crystal

@pytest.fixture
def square_2d():
    t1 = complexArray([0, 1])
    t2 = complexArray([1, 0])
    square_crystal = Crystal(t1, t2)
    yield square_crystal

@pytest.fixture
def rectangular_2d():
    t1 = complexArray([0, 1])
    t2 = complexArray([2, 0])
    square_crystal = Crystal(t1, t2)
    yield square_crystal

def test_raise_rank_error():
    t1 = complexArray([1])
    t2 = complexArray([2])
    with pytest.raises(ValueError) as err:
        crystal = Crystal(t1, t2, er=1, ur=1)

def test_1D_lattice_vectors(crystal_1D):
    T1_desired = 2 * pi * complexArray([1, 0])
    T1_calculated = crystal_1D.reciprocal_lattice_vectors[0]

    assert_almost_equal(T1_desired, T1_calculated);

def test_crystal_2D_more_dimensions():
    t1 = complexArray([1, 0, 0])
    t2 = complexArray([0, 1, 0])
    crystal = Crystal(t1, t2, er=1, ur=1)
    T1_desired = 2 * pi * complexArray([1, 0])
    T2_desired = 2 * pi * complexArray([0, 1])

    T1_calculated, T2_calculated = crystal.reciprocal_lattice_vectors
    assert_almost_equal(T1_desired, T1_calculated);
    assert_almost_equal(T2_desired, T2_calculated);

def testCalculateReciprocalLatticeVectors():
    # A simple cubic 2D lattice
    t1 = complexArray([1,0])
    t2 = complexArray([0,1])
    squareCrystal = Crystal(t1, t2)
    T1Actual = 2 * pi * complexArray([1,0])
    T2Actual = 2 * pi * complexArray([0,1])
    reciprocalLatticeVectorsActual = (T1Actual, T2Actual)
    reciprocalLatticeVectorsCalculated = squareCrystal.reciprocal_lattice_vectors

    assert_almost_equal(reciprocalLatticeVectorsActual, reciprocalLatticeVectorsCalculated)

    # A rectangular 2D lattice
    t1 = complexArray([2,0])
    t2 = complexArray([0,1])
    rectangularCrystal = Crystal(t1, t2)
    T1Actual = 1 * pi * complexArray([1 , 0])
    T2Actual = 2 * pi * complexArray([0 , 1])
    reciprocalLatticeVectorsActual = (T1Actual, T2Actual)
    reciprocalLatticeVectorsCalculated = rectangularCrystal.reciprocal_lattice_vectors

    assert_almost_equal(reciprocalLatticeVectorsActual, reciprocalLatticeVectorsCalculated);

def testDetermineCrystalType():
    # A square lattice
    t1 = complexArray([1,0]);
    t2 = complexArray([0,1]);
    squareCrystal = Crystal(t1, t2)
    crystalTypeActual = "SQUARE"
    crystalTypeCalculated = squareCrystal.crystalType
    assert_equal(crystalTypeActual, crystalTypeCalculated)


    # A rectangular lattice
    t1 = complexArray([1,0])
    t2 = complexArray([0,2])
    rectangularCrystal = Crystal(t1, t2)
    crystalTypeActual = "RECTANGULAR"
    crystalTypeCalculated = rectangularCrystal.crystalType
    assert_equal(crystalTypeActual, crystalTypeCalculated)

def test_lattice_vectors_3d(oblique_3d):
    desired_T1 = 2*np.pi * complexArray([1, 1, -3])
    desired_T2 = 2*np.pi * complexArray([0, 1, -1])
    desired_T3 = 2 * np.pi * complexArray([0, -1, 2])
    actual_T1 = oblique_3d.reciprocal_lattice_vectors[0]
    actual_T2 = oblique_3d.reciprocal_lattice_vectors[1]
    actual_T3 = oblique_3d.reciprocal_lattice_vectors[2]
    assert_allclose(desired_T1, actual_T1)
    assert_allclose(desired_T2, actual_T2)
    assert_allclose(desired_T3, actual_T3)

def test_impossible_crystal():
    with pytest.raises(ValueError):
        four_d_crystal = Crystal([1, 0, 0, 0],[0, 1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1])


def test_under_rank_crystal():
    with pytest.raises(ValueError):
        wrong_crystal = Crystal([0, 1], [2, 3], [4, 5])


def test_type_2d_oblique(oblique_2d):
    assert oblique_2d._crystal_type() == 'OBLIQUE'


def test_type_2d_square(square_2d):
    assert square_2d._crystal_type() == 'SQUARE'


def test_type_2d_rectangular(rectangular_2d):
    assert rectangular_2d._crystal_type() == 'RECTANGULAR'


def test_type_3d_square(square_3d):
    assert square_3d._crystal_type() == 'SQUARE'


def test_type_3d_rectangular(rectangular_3d):
    assert rectangular_3d._crystal_type() == 'RECTANGULAR'


def test_type_3d_oblique(oblique_3d):
    assert oblique_3d._crystal_type() == 'OBLIQUE'
