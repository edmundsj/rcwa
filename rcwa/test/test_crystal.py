from rcwa import Source, Crystal
from rcwa.shorthand import *
from rcwa.testing import *
from numpy.testing import assert_equal
import pytest

@pytest.fixture
def crystal_1D():
    t1 = complexArray([1, 0])
    square_crystal = Crystal(t1)
    return square_crystal

def test_raise_rank_error():
    t1 = complexArray([1])
    t2 = complexArray([2])
    with pytest.raises(ValueError) as err:
        crystal = Crystal(t1, t2, er=1, ur=1)

def test_1D_lattice_vectors(crystal_1D):
    T1_desired = 2 * pi * complexArray([1, 0])
    T1_calculated = crystal_1D.reciprocalLatticeVectors[0]

    assert_almost_equal(T1_desired, T1_calculated);

def test_crystal_2D_more_dimensions():
    t1 = complexArray([1, 0, 0])
    t2 = complexArray([0, 1, 0])
    crystal = Crystal(t1, t2, er=1, ur=1)
    T1_desired = 2 * pi * complexArray([1, 0])
    T2_desired = 2 * pi * complexArray([0, 1])

    T1_calculated, T2_calculated = crystal.reciprocalLatticeVectors
    assert_almost_equal(T1_desired, T1_calculated);
    assert_almost_equal(T2_desired, T2_calculated);

def test_1D_symmetry_points(crystal_1D):
    assert crystal_1D.symmetry_points is None
    assert crystal_1D.key_symmetry_names is None

def testCalculateReciprocalLatticeVectors():
    # A simple cubic 2D lattice
    t1 = complexArray([1,0])
    t2 = complexArray([0,1])
    squareCrystal = Crystal(t1, t2)
    T1Actual = 2 * pi * complexArray([1,0])
    T2Actual = 2 * pi * complexArray([0,1])
    reciprocalLatticeVectorsActual = (T1Actual, T2Actual)
    reciprocalLatticeVectorsCalculated = squareCrystal.reciprocalLatticeVectors

    assert_almost_equal(reciprocalLatticeVectorsActual, reciprocalLatticeVectorsCalculated)

    # A rectangular 2D lattice
    t1 = complexArray([2,0])
    t2 = complexArray([0,1])
    rectangularCrystal = Crystal(t1, t2)
    T1Actual = 1 * pi * complexArray([1 , 0])
    T2Actual = 2 * pi * complexArray([0 , 1])
    reciprocalLatticeVectorsActual = (T1Actual, T2Actual)
    reciprocalLatticeVectorsCalculated = rectangularCrystal.reciprocalLatticeVectors

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

def testGenerateKeySymmetryPoints():

    # A square lattice
    t1 = complexArray([1,0])
    t2 = complexArray([0,1])
    squareCrystal = Crystal(t1, t2)
    T1 = 2*pi*complexArray([1, 0])
    T2 = 2*pi*complexArray([0,1])

    keySymmetryPointsActual = [0.5 * T1, 0*T1, 0.5 * (T1 + T2)]
    keySymmetryNamesActual = ["X", "G", "M"]
    keySymmetryPointsCalculated = squareCrystal.symmetry_points
    keySymmetryNamesCalculated = squareCrystal.key_symmetry_names

    assert_equal(keySymmetryPointsActual, keySymmetryPointsCalculated);
    assert_equal(keySymmetryNamesActual, keySymmetryNamesCalculated);

    # A rectangular Lattice
    t1 = complexArray([1,0])
    t2 = complexArray([0,2])
    rectangularCrystal = Crystal(t1, t2)
    T1 = 2*pi*complexArray([1, 0]);
    T2 = pi*complexArray([0,1]);

    keySymmetryPointsActual = [0.5 * T1, 0 * T1, 0.5 * T2, 0.5 * (T1 + T2)];
    keySymmetryNamesActual = ["X", "G", "Y", "S"];
    keySymmetryPointsCalculated = rectangularCrystal.symmetry_points;
    keySymmetryNamesCalculated = rectangularCrystal.key_symmetry_names;

    assert_equal(keySymmetryPointsActual, keySymmetryPointsCalculated);
    assert_equal(keySymmetryNamesActual, keySymmetryNamesCalculated);
