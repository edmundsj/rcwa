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

from matrices import *
from fresnel import *
from convolution import generateConvolutionMatrix
from shorthandTest import *
from matrixParser import *

class Test:
    def __init__(self):
        self.messages = []; # Messages sent back from our tests (strings)
        self.statuses = []; # statuses sent back from our tests (boolean)
        self.unitTestsEnabled = True;
        self.integrationTestsEnabled = True;

    def printResults(self):
        for s, i in zip(self.statuses, range(len(self.statuses))):
            if(s == False):
                print(self.messages[i]);
        print(f"{self.statuses.count(True)} PASSED, {self.statuses.count(False)} FAILED");

    def testCaller(self, testFunction, *args):
        test_status = False; # By default assume we failed the test.
        test_message = f"{testFunction.__name__}({args}): ";

        try:
            print(f"Calling function {testFunction.__name__} ... ", end=" ");
            testFunction(*args);
            print("OK");
            test_status = True;
            self.statuses.append(test_status);
            self.messages.append(test_message);
        except AssertionError as ae:
            print("FAIL");
            test_message += "FAILED";
            test_message += str(ae);
            self.statuses.append(test_status);
            self.messages.append(test_message);

    def runUnitTests(self):
        print("--------- RUNNING UNIT TESTS... ----------");
        self.testCaller(self.testSetupData1x1Harmonics)
        self.testCaller(self.testGenerateConvolutionMatrix)
        print("--------- END UNIT TESTS... ----------");

    def runIntegrationTests(self):
        """
        Runs integration tests to verify s-parameters for composite code, to verify the output field
        for a given input field, and to verify the reflectance/transmittance and enforce power 
        conservation.
        """

        print("--------- RUNNING INTEGRATION TESTS... ----------");
        self.testCaller(self.itestGlobalScatteringMatrix);

        print("--------- END INTEGRATION TESTS... ----------");

    def testSetupData1x1Harmonics(self):
        self.setupData1x1Harmonics()
        absoluteTolerance = 1e-4
        relativeTolerance = 1e-3
        lambda0Actual = 0.02
        thetaActual = 0
        phiActual = 0

        pTEActual = 1
        pTMActual = 0

        urReflectionRegionActual = 1.0
        erReflectionRegionActual = 2.0
        urTransmissionRegionActual = 1.0
        erTransmissionRegionActual = 9.0
        urDeviceRegionActual = 1.0
        erDeviceRegionActual = 6.0

        xPeriodActual = 0.0175
        yPeriodActual = 0.015
        layer1ThicknessActual = 0.005
        layer2ThicknessActual = 0.003
        triangleWidthActual = 0.012

        NxActual = 512
        NyActual = 439
        xHarmonicsActual = 1
        yHarmonicsActual = 1
        dxActual = 3.418e-5
        dyActual = 3.4169e-5

        assertEqual(lambda0Actual, self.lambda0, "Free space wavelength not equal")
        assertEqual(thetaActual, self.theta, "Angle theta not equal")
        assertEqual(phiActual, self.phi, "Angle phi not equal")
        assertEqual(pTEActual, self.pTE, "TE polarization amount not equal")
        assertEqual(pTMActual, self.pTM, "TM polarization amount not equal")
        assertEqual(urReflectionRegionActual, self.urReflectionRegion, "ur in reflection region not equal")
        assertEqual(erReflectionRegionActual, self.erReflectionRegion, "er in reflection region not equal")
        assertEqual(erTransmissionRegionActual, self.erTransmissionRegion,
                "er in transmission region not equal")
        assertEqual(urTransmissionRegionActual, self.urTransmissionRegion,
                "ur in transmission region not equal")
        assertEqual(erDeviceRegionActual, self.erDeviceRegion, "er in device region not equal")
        assertEqual(urDeviceRegionActual, self.urDeviceRegion, "ur in device region not equal")
        assertAlmostEqual(NxActual, self.Nx, absoluteTolerance, relativeTolerance,
                "Nx not equal")
        assertAlmostEqual(NyActual, self.Ny, absoluteTolerance, relativeTolerance,
                "Ny not equal")
        assertAlmostEqual(xHarmonicsActual, self.numberHarmonics[0],
                errorMessage="number x harmonics not equal")
        assertEqual(yHarmonicsActual, self.numberHarmonics[1],
                errorMessage="number y harmonics not equal")
        assertAlmostEqual(xPeriodActual, self.xPeriod, absoluteTolerance, relativeTolerance,
                "x Period not equal")
        assertAlmostEqual(yPeriodActual, self.yPeriod, absoluteTolerance, relativeTolerance,
                "y Period not equal")
        assertAlmostEqual(layer1ThicknessActual, self.layerThickness[0], absoluteTolerance, relativeTolerance,
                "Layer 1 thicknesses not equal")
        assertAlmostEqual(layer2ThicknessActual, self.layerThickness[1], absoluteTolerance, relativeTolerance,
                "Layer 2 thickness not equal")
        assertAlmostEqual(dxActual, self.dx, absoluteTolerance, relativeTolerance,
                "dx not equal")
        assertAlmostEqual(dyActual, self.dy, absoluteTolerance, relativeTolerance,
                "dy not equal")

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

    def itestReflectanceTransmittance(self):
        """
        Tests that the reflectance and transmittance are correct, and that conservation of power
        is enforced.
        """

        RActual = 0.4403;
        TActual = 0.5597;
        CONActual = 1;

    def setupData1x1Harmonics(self):
        cm = 1e-2
        deg = pi / 180
        self.lambda0 = 2 * cm
        self.theta = 0 * deg
        self.phi = 0 * deg
        self.pTE = 1
        self.pTM = 0

        self.urReflectionRegion = 1.0
        self.erReflectionRegion = 2.0
        self.urTransmissionRegion = 1.0
        self.erTransmissionRegion = 9.0
        self.erDeviceRegion = 6.0
        self.urDeviceRegion = 1.0
        self.xPeriod = 1.75 * cm
        self.yPeriod = 1.5 * cm
        thicknessLayer1 = 0.5 * cm
        thicknessLayer2 = 0.3 * cm

        self.Nx = 512;
        self.Ny = round(self.Nx * self.yPeriod / self.xPeriod);
        (spatialHarmonicsX, spatialHarmonicsY) = (1, 1)
        self.numberHarmonics = (spatialHarmonicsX, spatialHarmonicsY)

        self.dx = self.xPeriod / self.Nx;
        self.dy = self.yPeriod / self.Ny;
        self.xCoors = np.linspace(-self.xPeriod/2 + self.dx/2, self.xPeriod/2 - self.dx/2, self.Nx)
        self.yCoors = np.linspace(-self.yPeriod/2 + self.dy/2, self.yPeriod/2 - self.dy/2, self.Ny)

        self.triangleWidth = 0.8 * self.xPeriod
        self.triangleHeight = 0.5 * sqrt(3) * self.triangleWidth

        # NOTE - this assumes that our matrices have the x component as the 2nd index and the y component
        # as the third, for ease of indexing. [layer][x][y]
        self.urData = self.urDeviceRegion * complexOnes((2, self.Nx, self.Ny))
        self.erData = self.erDeviceRegion * complexOnes((2, self.Nx, self.Ny))
        triangleData = np.transpose(np.loadtxt('test/triangleData.csv', delimiter=','))
        self.erData[0] = triangleData;
        self.layerThickness = [thicknessLayer1, thicknessLayer2]

    test_class = Test(); # Create a new test class
    if(test_class.unitTestsEnabled == True):
        test_class.runUnitTests();
    if(test_class.integrationTestsEnabled == True):
        test_class.runIntegrationTests();
    test_class.printResults();

main();
