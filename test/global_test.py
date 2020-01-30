# TODO:
# 0. Clean up your code by doing:
#
# 1. Write unit tests for the following methods:
#
# 2. Write integration tests for the following:

import sys
sys.path.append('core');
sys.path.append('test')
from matplotlib import pyplot as plt

import unittest
from matrices import *
from fresnel import *
from convolution import generateConvolutionMatrix
from shorthandTest import *
from matrixParser import *

class Test(unittest.TestCase):
    def testGenerateConvolutionMatrix(self):
        absoluteTolerance = 1e-4
        relativeTolerance = 1e-3

        self.setupData1x1Harmonics()
        convolutionMatrixCalculated = generateConvolutionMatrix(self.urData[0], self.numberHarmonics)
        convolutionMatrixActual = 1
        assertAlmostEqual(convolutionMatrixActual, convolutionMatrixCalculated, absoluteTolerance,
                relativeTolerance, "UR convolution matrices for layer 1 not equal")


        convolutionMatrixCalculated = generateConvolutionMatrix(self.erData[0], self.numberHarmonics)
        convolutionMatrixActual = 5.0449
        assertAlmostEqual(convolutionMatrixActual, convolutionMatrixCalculated, absoluteTolerance,
                relativeTolerance, "ER convolution matrices for layer 1 not equal")

        convolutionMatrixActual = 1
        convolutionMatrixCalculated = generateConvolutionMatrix(self.urData[1], self.numberHarmonics)
        assertAlmostEqual(convolutionMatrixActual, convolutionMatrixCalculated, absoluteTolerance,
                relativeTolerance, "UR convolution matrices for layer 2 not equal")

        convolutionMatrixActual = 6
        convolutionMatrixCalculated = generateConvolutionMatrix(self.erData[1], self.numberHarmonics)
        assertAlmostEqual(convolutionMatrixActual, convolutionMatrixCalculated, absoluteTolerance,
                relativeTolerance, "ER convolution matrices for layer 2 not equal")

    def itestGlobalScatteringMatrix(self):
        """
        Tests that the global scattering matrix is correct for the overall structure. If this works,
        then pretty much everything should work.
        """
        absoluteTolerance = 0.0001;
        relativeTolerance = 0.001;
        pTE = 1/sqrt(2);
        pTM = (1j)/sqrt(2);

        l0 = 2.7;
        k0 = 2*np.pi / l0;
        thetaInDegrees = 57;
        phiInDegrees = 23;
        theta = np.pi / 180.0 * thetaInDegrees;
        phi = np.pi / 180.0 * phiInDegrees;

        er = [2.0, 1.0];
        ur = [1.0, 3.0];
        L = [0.25*l0, 0.5*l0];

        erReflectionRegion = 1.4;
        urReflectionRegion = 1.2;
        erTransmissionRegion = 1.8;
        urTransmissionRegion = 1.6;

        # First, calculate the incident k-vector
        kVector = calculateKVector(theta, phi, erReflectionRegion, urReflectionRegion);
        kx = kVector[0];
        ky = kVector[1];
        kzReflectionRegion = kVector[2];

        # Calculate gap medium parameters
        erGap = 1 + sq(kx) + sq(ky); # This can be anything, but this simplifies an intermediate matrix
        urGap = 1;
        kzGap = calculateKz(kx, ky, erGap, urGap); # Should be 1.
        (Vg, Wg) = calculateVWXMatrices(kx, ky, kzGap, erGap, urGap);
        # THIS PART LOOKS GOOD.

        # Initialize the global scattering matrix
        Sglobal = complexZeros((2,2,2,2));
        Sglobal[1,0] = complexIdentity(2);
        Sglobal[0,1] = complexIdentity(2);

        # Now, loop through the layers - THIS PART MATCHES WHAT WE WANT.
        # BOTH THE SGLOBAL MATRICES AND THE SI MATRICES MATCH WHAT THEY SHOULD.
        numberOfLayers = len(L);
        for i in range(numberOfLayers):
            Si = calculateInternalSMatrix(kx, ky, er[i], ur[i], k0, L[i], Wg, Vg);

            # Finally, compute the redheffer product with the current global matrix
            Sglobal = calculateRedhefferProduct(Sglobal, Si);

        # Calculate the reflection and transmission region s matrices
        SReflectionRegion = calculateReflectionRegionSMatrix(kx, ky,
                erReflectionRegion, urReflectionRegion, Wg, Vg);
        # THE TRANSMISSION REGION MATRIX LOOKS WRONG.
        STransmissionRegion = calculateTransmissionRegionSMatrix(kx, ky,
                erTransmissionRegion, urTransmissionRegion, Wg, Vg);

        # Finally, compute the redheffer star product to connect our global matrix to the external
        # regions
        Sglobal = calculateRedhefferProduct(Sglobal, STransmissionRegion);
        Sglobal = calculateRedhefferProduct(SReflectionRegion, Sglobal);

        SGlobalCalculated = Sglobal;

        SGlobalActual = complexZeros((2,2,2,2));
        SGlobalActual[0,0] = complexArray([
            [-0.6018 + 0.3062j, -0.0043 + 0.0199j],
            [-0.0043 + 0.0199j, -0.5935 + 0.2678j]]);
        SGlobalActual[0,1] = complexArray([
            [0.5766 - 0.3110j, -0.0919 + 0.0469j],
            [-0.0919 + 0.0469j, 0.7542 - 0.4016j]]);
        SGlobalActual[1,0] = complexArray([
            [0.7415 - 0.4007j, 0.0716 - 0.0409j],
            [0.0716 - 0.0409j, 0.6033 - 0.3218j]]);
        SGlobalActual[1,1] = complexArray([
            [0.5861 - 0.3354j, 0.0170 + 0.0042j],
            [0.0170 + 0.0042j, 0.5533 - 0.3434j]]);

        assertAlmostEqual(SGlobalCalculated, SGlobalActual, absoluteTolerance, relativeTolerance);

if __name__ == '__main__':
    unittest.main()
