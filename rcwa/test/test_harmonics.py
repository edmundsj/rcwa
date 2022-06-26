import unittest
from rcwa import Crystal, Source
from rcwa.shorthand import *
from rcwa.testing import *
from rcwa.harmonics import *
import numpy as np

class testHarmonicFunctions(unittest.TestCase):
    def testGetXComponents(self):
        testVector1 = complexArray([0.163, 0.5, 0.888]);
        testVector2 = complexArray([0.246, 0.99, 0.2]);
        xComponentsCalculated = x_components(testVector1, testVector2);
        xComponentsActual = [0.163, 0.246];
        assert_almost_equal(xComponentsActual, xComponentsCalculated);

        testVector1 = complexArray([[0.183], [0.5], [0.888]]);
        testVector2 = complexArray([[0.266], [0.99], [0.2]]);
        xComponentsCalculated = x_components(testVector1, testVector2);
        xComponentsActual = [0.183, 0.266];
        assert_almost_equal(xComponentsActual, xComponentsCalculated);

        testVector1 = complexArray([[1.173], [0.7]]);
        testVector2 = complexArray([1.256, 1.99]);
        xComponentsCalculated = x_components(testVector1, testVector2);
        xComponentsActual = [1.173, 1.256];
        assert_almost_equal(xComponentsActual, xComponentsCalculated);

    def testGetYComponents(self):
        testVector1 = complexArray([0.173, 0.4, 0.888]);
        testVector2 = complexArray([0.256, 0.89, 0.2]);
        yComponentsCalculated = y_components(testVector1, testVector2);
        yComponentsActual = [0.4, 0.89];
        assert_almost_equal(yComponentsActual, yComponentsCalculated);

        testVector1 = complexArray([[0.173], [0.5], [0.888]]);
        testVector2 = complexArray([[0.256], [0.99], [0.2]]);
        yComponentsCalculated = y_components(testVector1, testVector2);
        yComponentsActual = [0.5, 0.99];
        assert_almost_equal(yComponentsActual, yComponentsCalculated);

        testVector1 = complexArray([[0.173], [0.7]]);
        testVector2 = complexArray([0.256, 1.99]);
        yComponentsCalculated = y_components(testVector1, testVector2);
        yComponentsActual = [0.7, 1.99];
        assert_almost_equal(yComponentsActual, yComponentsCalculated);

    def testCalculateZeroHarmonicLocation(self):
        harmonicNumber1 = 5;
        harmonicNumber2 = 6;
        numberHarmonics = (harmonicNumber1, harmonicNumber2)
        zeroHarmonicLocationsCalculated = zero_harmonic(numberHarmonics)
        zeroHarmonicLocationsActual = [2, 3];
        assert_almost_equal(zeroHarmonicLocationsActual, zeroHarmonicLocationsCalculated);

    def testCalculateMinHarmonic(self):
        harmonicNumber1 = 5;
        harmonicNumber2 = 6;
        numberHarmonics = (harmonicNumber1, harmonicNumber2)
        minHarmonicCalculated= min_harmonic(numberHarmonics)
        minHarmonicActual = [-2, -3];
        assert_almost_equal(minHarmonicActual, minHarmonicCalculated);

    def testCalculateMaxHarmonic(self):
        harmonicNumber1 = 5;
        harmonicNumber2 = 6;
        numberHarmonics = (harmonicNumber1, harmonicNumber2)
        maxHarmonicCalculated= max_harmonic(numberHarmonics)
        maxHarmonicActual = [2, 2];
        assert_almost_equal(maxHarmonicActual, maxHarmonicCalculated);

    # ALL THESE TESTS WILL FAIL BECAUSE I NOW NEED TO PASS IN A SOURCE WHICH HAS AN UNKNOWN K0.
    def testGenerateKxMatrix(self):
        absoluteTolerance = 1e-4;
        relativeTolerance = 1e-3;

        # Test our KX matrix at the gamma point
        kxMatrixActual = self.KxMatrixGPoint;
        sourceGPoint = Source()
        sourceGPoint._k_incident = self.incidentKVectorGPoint
        sourceGPoint.k0 = 1
        kxMatrixCalculated = kx_matrix(sourceGPoint, self.crystal, self.numberHarmonics)
        assert_almost_equal(kxMatrixActual, kxMatrixCalculated, absoluteTolerance, relativeTolerance,
                "testHarmonics: Kx at Gamma Point");

        # Test our KX matrix at the X point
        kxMatrixActual = self.KxMatrixXPoint;
        sourceXPoint = Source()
        sourceXPoint._k_incident = self.incidentKVectorXPoint
        sourceXPoint.k0 = 1
        kxMatrixCalculated = kx_matrix(sourceXPoint, self.crystal, self.numberHarmonics)
        assert_almost_equal(kxMatrixActual, kxMatrixCalculated, absoluteTolerance, relativeTolerance,
                "testHarmonics: Kx at X Point");

        # Test our KX matrix at the M point
        kxMatrixActual = self.KxMatrixMPoint;
        sourceMPoint = Source()
        sourceMPoint._k_incident = self.incidentKVectorMPoint
        sourceMPoint.k0 = 1
        kxMatrixCalculated = kx_matrix(sourceMPoint, self.crystal, self.numberHarmonics)
        assert_almost_equal(kxMatrixActual, kxMatrixCalculated, absoluteTolerance, relativeTolerance,
                "testHarmonics: Kx at M Point");

    def testGenerateKyMatrix(self):
        absoluteTolerance = 1e-4;
        relativeTolerance = 1e-3;

        # Test our KY matrix at the gamma point
        kyMatrixActual = self.KyMatrixGPoint;
        sourceGPoint = Source()
        sourceGPoint._k_incident = self.incidentKVectorGPoint
        sourceGPoint.k0 = 1
        kyMatrixCalculated = ky_matrix(sourceGPoint, self.crystal, self.numberHarmonics)
        assert_almost_equal(kyMatrixActual, kyMatrixCalculated, absoluteTolerance, relativeTolerance,
                "testHarmonics: Ky at Gamma Point");

        # Test our KY matrix at the X point
        kyMatrixActual = self.KyMatrixXPoint;
        sourceXPoint = Source()
        sourceXPoint._k_incident = self.incidentKVectorXPoint
        sourceXPoint.k0 = 1
        kyMatrixCalculated = ky_matrix(sourceXPoint, self.crystal, self.numberHarmonics)
        assert_almost_equal(kyMatrixActual, kyMatrixCalculated, absoluteTolerance, relativeTolerance,
                "testHarmonics: Ky at X Point");

        # Test our KY matrix at the M point
        kyMatrixActual = self.KyMatrixMPoint;
        sourceMPoint = Source()
        sourceMPoint._k_incident = self.incidentKVectorMPoint
        sourceMPoint.k0 = 1
        kyMatrixCalculated = ky_matrix(sourceMPoint, self.crystal, self.numberHarmonics)
        assert_almost_equal(kyMatrixActual, kyMatrixCalculated, absoluteTolerance, relativeTolerance,
                "testHarmonics: Ky at M Point");

    def setUp(self):
        self.numberHarmonics = (3, 3, 1)
        self.matrixDimensions = np.prod(self.numberHarmonics);
        self.matrixShape = (self.matrixDimensions, self.matrixDimensions);

        a = 1;
        self.ax = a;
        self.ay = a;
        self.r = 0.35 * a;
        self.er = 9.0;

        self.Nx = 512;
        self.Ny = 512;
        self.dx = self.ax / self.Nx;
        self.dy = self.ay / self.Ny;

        self.t1 = complexArray([self.ax, 0])
        self.t2 = complexArray([0, self.ay])
        self.T1 = (2 * pi / self.ax) * complexArray([1, 0]);
        self.T2 = (2 * pi / self.ay) * complexArray([0, 1]);

        self.incidentKVectorGPoint= 0*self.T1;
        self.incidentKVectorXPoint = 0.5*self.T1;
        self.incidentKVectorMPoint = 0.5*self.T1 + 0.5*self.T2;

        xcoors = np.linspace(-self.ax/2 + self.dx/2, self.ax/2 - self.dx/2, self.Nx);
        ycoors = np.linspace(-self.ay/2 + self.dy/2, self.ay/2 - self.dy/2, self.Ny);
        (X, Y) = np.meshgrid(xcoors, ycoors);
        self.UR = complexOnes((self.Nx, self.Ny));
        self.ER = (self.er-1) * np.heaviside(sq(X) + sq(Y) - sq(self.r),1)
        self.ER = self.ER + 1;
        source = Source()
        self.crystal = Crystal(self.t1, self.t2, er=self.ER, ur=self.UR)

        # The data for Kx, Ky, and Kz will be re-used at each point of key symmetry
        self.KxMatrixGPoint = complexZeros(self.matrixShape);
        self.KxMatrixGPoint[0,0] = 6.2832;
        self.KxMatrixGPoint[2,2] = -6.2832;
        self.KxMatrixGPoint[3,3] = 6.2832;
        self.KxMatrixGPoint[5,5] = -6.2832;
        self.KxMatrixGPoint[6,6] = 6.2832;
        self.KxMatrixGPoint[8,8] = -6.2832;

        self.KyMatrixGPoint = complexZeros(self.matrixShape);
        self.KyMatrixGPoint[0,0] = 6.2832;
        self.KyMatrixGPoint[1,1] = 6.2832;
        self.KyMatrixGPoint[2,2] = 6.2832;
        self.KyMatrixGPoint[6,6] = -6.2832;
        self.KyMatrixGPoint[7,7] = -6.2832;
        self.KyMatrixGPoint[8,8] = -6.2832;

        self.KxMatrixXPoint = complexZeros(self.matrixShape);
        self.KxMatrixXPoint[0,0] = 9.4248;
        self.KxMatrixXPoint[1,1] = 3.1416;
        self.KxMatrixXPoint[2,2] = -3.1416;
        self.KxMatrixXPoint[3,3] = 9.4248;
        self.KxMatrixXPoint[4,4] = 3.1416;
        self.KxMatrixXPoint[5,5] = -3.1416;
        self.KxMatrixXPoint[6,6] = 9.4248;
        self.KxMatrixXPoint[7,7] = 3.1416;
        self.KxMatrixXPoint[8,8] = -3.1416;

        diagonalValuesKyXPoint = complexArray([6.2832, 6.2832, 6.2832, 0, 0, 0, -6.2832, -6.2832, -6.2832]);
        self.KyMatrixXPoint = np.diag(diagonalValuesKyXPoint);

        diagonalValuesKxMPoint = complexArray([9.4248, 3.1416, -3.1416, 9.4248, 3.1416,
            -3.1416, 9.4248, 3.1416, -3.1416]);
        self.KxMatrixMPoint = np.diag(diagonalValuesKxMPoint);

        diagonalValuesKyMPoint = complexArray([9.4248, 9.4248, 9.4248, 3.1416, 3.1416,
            3.1416,-3.1416, -3.1416, -3.1416]);
        self.KyMatrixMPoint = np.diag(diagonalValuesKyMPoint);

if __name__ == 'main':
    unittest.main()
