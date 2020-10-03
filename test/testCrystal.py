import context
from source import *
from crystal import *
import unittest
from shorthandTest import *

class Test(unittest.TestCase):

    def testCalculateReciprocalLatticeVectors(self):
        # A simple cubic 2D lattice
        t1 = complexArray([1,0]);
        t2 = complexArray([0,1]);
        squareCrystal = Crystal(1, 1, t1, t2)
        T1Actual = 2 * pi * complexArray([1,0]);
        T2Actual = 2 * pi * complexArray([0,1]);
        reciprocalLatticeVectorsActual = (T1Actual, T2Actual);
        reciprocalLatticeVectorsCalculated = squareCrystal.reciprocalLatticeVectors

        assertAlmostEqual(reciprocalLatticeVectorsActual, reciprocalLatticeVectorsCalculated);

        # A rectangular 2D lattice
        t1 = complexArray([2,0]);
        t2 = complexArray([0,1]);
        rectangularCrystal = Crystal(1, 1, t1, t2);
        T1Actual = 1 * pi * complexArray([1 , 0]);
        T2Actual = 2 * pi * complexArray([0 , 1]);
        reciprocalLatticeVectorsActual = (T1Actual, T2Actual);
        reciprocalLatticeVectorsCalculated = rectangularCrystal.reciprocalLatticeVectors

        assertAlmostEqual(reciprocalLatticeVectorsActual, reciprocalLatticeVectorsCalculated);

    def testDetermineCrystalType(self):
        # A square lattice
        t1 = complexArray([1,0]);
        t2 = complexArray([0,1]);
        squareCrystal = Crystal(1, 1, t1, t2)
        crystalTypeActual = "SQUARE"
        crystalTypeCalculated = squareCrystal.crystalType
        assertStringEqual(crystalTypeActual, crystalTypeCalculated)


        # A rectangular lattice
        t1 = complexArray([1,0]);
        t2 = complexArray([0,2]);
        rectangularCrystal = Crystal(1, 1, t1, t2)
        crystalTypeActual = "RECTANGULAR";
        crystalTypeCalculated = rectangularCrystal.crystalType
        assertStringEqual(crystalTypeActual, crystalTypeCalculated)

    def testGenerateKeySymmetryPoints(self):

        # A square lattice
        t1 = complexArray([1,0]);
        t2 = complexArray([0,1]);
        squareCrystal = Crystal(1, 1, t1, t2)
        T1 = 2*pi*complexArray([1, 0]);
        T2 = 2*pi*complexArray([0,1]);

        keySymmetryPointsActual = [0.5 * T1, 0*T1, 0.5 * (T1 + T2)]
        keySymmetryNamesActual = ["X", "G", "M"]
        keySymmetryPointsCalculated = squareCrystal.keySymmetryPoints
        keySymmetryNamesCalculated = squareCrystal.keySymmetryNames

        assertArrayEqual(keySymmetryPointsActual, keySymmetryPointsCalculated);
        assertArrayEqual(keySymmetryNamesActual, keySymmetryNamesCalculated);

        # A rectangular Lattice
        t1 = complexArray([1,0])
        t2 = complexArray([0,2])
        rectangularCrystal = Crystal(1, 1, t1, t2)
        T1 = 2*pi*complexArray([1, 0]);
        T2 = pi*complexArray([0,1]);

        keySymmetryPointsActual = [0.5 * T1, 0 * T1, 0.5 * T2, 0.5 * (T1 + T2)];
        keySymmetryNamesActual = ["X", "G", "Y", "S"];
        keySymmetryPointsCalculated = rectangularCrystal.keySymmetryPoints;
        keySymmetryNamesCalculated = rectangularCrystal.keySymmetryNames;

        assertArrayEqual(keySymmetryPointsActual, keySymmetryPointsCalculated);
        assertArrayEqual(keySymmetryNamesActual, keySymmetryNamesCalculated);

if __name__ == '__main__':
    unittest.main()
