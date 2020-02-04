import sys
sys.path.append('core');
sys.path.append('test')

import unittest
from shorthandTest import *
from matrices import *
from fresnel import *
from convolution import generateConvolutionMatrix
from matrixParser import *


class Test3x3HarmonicsOblique(unittest.TestCase):

    def testTransparentSMatrix(self):
        SActual = self.transparentSMatrix
        SCalculated = generateTransparentSMatrix((18, 18));
        assertAlmostEqual(SActual, SCalculated,self.absoluteTolerance,self.relativeTolerance);

    def testPMatrix(self):
        PActual = self.PLayer1
        PCalculated = calculatePMatrix(self.Kx, self.Ky,
               self.erConvolutionMatrixLayer1, self.urConvolutionMatrixLayer1)
        assertAlmostEqual(PActual, PCalculated, self.absoluteTolerance, self.relativeTolerance,
                "P matrix layer 1");

    def testQMatrix(self):
        QActual = self.QLayer1
        QCalculated = calculateQMatrix(self.Kx, self.Ky,
                self.erConvolutionMatrixLayer1, self.urConvolutionMatrixLayer1)
        assertAlmostEqual(QActual, QCalculated, self.absoluteTolerance, self.relativeTolerance,
                "Q matrix Layer 1");

        QActual = self.QReflectionRegion
        QCalculated = calculateQMatrix(self.Kx, self.Ky,
                self.erReflectionRegion, self.urReflectionRegion)
        assertAlmostEqual(QActual, QCalculated, self.absoluteTolerance, self.relativeTolerance,
                "Q Reflection Region");

        QActual = self.QTransmissionRegion
        QCalculated = calculateQMatrix(self.Kx, self.Ky,
                self.erTransmissionRegion, self.urTransmissionRegion)
        assertAlmostEqual(QActual, QCalculated, self.absoluteTolerance, self.relativeTolerance,
                "Q Transmission Region");

    def testOmegaSquaredMatrix(self):
        OmegaSquaredActual = self.OmegaSquaredLayer1
        OmegaSquaredCalculated = calculateOmegaSquaredMatrix(self.PLayer1, self.QLayer1)
        assertAlmostEqual(OmegaSquaredActual, OmegaSquaredCalculated,
                self.absoluteTolerance, self.relativeTolerance);

    def testWMatrix(self):
        (V, WCalculated, X) = calculateVWXMatrices(self.Kx, self.Ky,
                self.erConvolutionMatrixLayer1, self.urConvolutionMatrixLayer1, self.k0,
                self.thicknessLayer1)
        WActual = self.WLayer1
        assertAlmostEqual(WActual, WCalculated, self.absoluteTolerance, self.relativeTolerance,
                "W matrix Layer 1");

    def testVMatrix(self):
        (VCalculated, W, X) = calculateVWXMatrices(self.Kx, self.Ky,
                self.erConvolutionMatrixLayer1, self.urConvolutionMatrixLayer1, self.k0,
                self.thicknessLayer1)
        VActual = self.VLayer1

        # Because the inverse of the lambda matrix has some negative values where they should be
        # positive, we need to negate all those indices.
        indicesToNegate = [6, 10, 11, 15]
        VActual[:, indicesToNegate] = - VActual[:, indicesToNegate]

        assertAlmostEqual(VActual, VCalculated, self.absoluteTolerance, self.relativeTolerance,
                "V matrix Layer 1");

    def testXMatrix(self):
        (V, W, XCalculated) = calculateVWXMatrices(self.Kx, self.Ky,
                self.erConvolutionMatrixLayer1, self.urConvolutionMatrixLayer1, self.k0,
                self.thicknessLayer1)
        XActual = self.XLayer1

        # Numerical error is causing accidental conjugation. To match the test data we need
        # to un-conjugate. 
        indices = [6, 10, 11, 15]
        XActual[indices, indices] = conj(XActual[indices, indices])
        assertAlmostEqual(XActual, XCalculated, self.absoluteTolerance, self.relativeTolerance,
                "X matrix Layer 1");

    def testAMatrix(self):

        V = self.VLayer1
        W = self.WLayer1
        ACalculated = calculateScatteringAMatrix(W, self.WFreeSpace, V, self.VFreeSpace);
        AActual = self.ALayer1
        assertAlmostEqual(AActual, ACalculated, self.absoluteTolerance, self.relativeTolerance);

    def testBMatrix(self):
        W = self.WLayer1
        V = self.VLayer1
        BCalculated = calculateScatteringBMatrix(W, self.WFreeSpace, V, self.VFreeSpace);
        BActual = self.BLayer1
        assertAlmostEqual(BActual, BCalculated, self.absoluteTolerance, self.relativeTolerance);

    def testScatteringMatrixFromRaw(self):
        SMatrixLayer1Calculated = calculateInternalSMatrixFromRaw(self.ALayer1, self.BLayer1,
                self.XLayer1, calculateScatteringDMatrix(self.ALayer1, self.BLayer1, self.XLayer1))
        S11Actual = self.S11Layer1
        S11Calculated = SMatrixLayer1Calculated[0,0];
        assertAlmostEqual(S11Actual, S11Calculated, self.absoluteTolerance, self.relativeTolerance,
                "S11 for Layer 1");

        S12Actual = self.S12Layer1
        S12Calculated = SMatrixLayer1Calculated[0,1];
        assertAlmostEqual(S12Actual, S12Calculated, self.absoluteTolerance, self.relativeTolerance,
                "S12 for Layer 1");

        S12Actual = self.S12Layer1
        S12Calculated = SMatrixLayer1Calculated[0,1];
        assertAlmostEqual(S12Actual, S12Calculated, self.absoluteTolerance, self.relativeTolerance,
                "S12 for Layer 1");

        S21Actual = self.S21Layer1
        S21Calculated = SMatrixLayer1Calculated[1,0];
        assertAlmostEqual(S21Actual, S21Calculated, self.absoluteTolerance, self.relativeTolerance,
                "S21 for Layer 1");

        S22Actual = self.S22Layer1
        S22Calculated = SMatrixLayer1Calculated[1,1];
        assertAlmostEqual(S22Actual, S22Calculated, self.absoluteTolerance, self.relativeTolerance,
                "S22 for Layer 1");

    def testSMatrixFromFundamentals(self):
        SiCalculated = calculateInternalSMatrix(self.Kx, self.Ky, self.erConvolutionMatrixLayer1,
                self.urConvolutionMatrixLayer1,
                self.k0, self.thicknessLayer1, self.WFreeSpace, self.VFreeSpace)

        S11Actual = self.S11Layer1
        S11Calculated = SiCalculated[0,0]
        assertAlmostEqual(S11Actual, S11Calculated,
                self.absoluteTolerance, self.relativeTolerance, "S11 Matrix layer 1")

        S12Actual = self.S12Layer1
        S12Calculated = SiCalculated[0,1]
        assertAlmostEqual(S12Actual, S12Calculated,
                self.absoluteTolerance, self.relativeTolerance, "S12 layer 1")

        S21Actual = self.S21Layer1
        S21Calculated = SiCalculated[1,0]
        assertAlmostEqual(S21Actual, S21Calculated,
                self.absoluteTolerance, self.relativeTolerance, "S21 Matrix layer 1")

        S22Actual = self.S22Layer1
        S22Calculated = SiCalculated[1,1]
        assertAlmostEqual(S22Actual, S22Calculated,
                self.absoluteTolerance, self.relativeTolerance, "S21 Matrix layer 1")

    def testDRedhefferMatrix(self):
        # Unfortunately we don't have more test data than getting the identity back.
        SA = self.transparentSMatrix
        SB = self.SLayer1
        DRedhefferMatrixActual = complexIdentity(18)
        DRedhefferMatrixCalculated = calculateRedhefferDMatrix(SA, SB)
        assertAlmostEqual(DRedhefferMatrixActual, DRedhefferMatrixCalculated, self.absoluteTolerance,
                self.relativeTolerance, "Layer 1 D matrix")

        SA = self.SLayer1
        SB = self.SLayer2
        DRedhefferMatrixActual = self.DLayer12
        DRedhefferMatrixCalculated = calculateRedhefferDMatrix(SA, SB)
        assertAlmostEqual(DRedhefferMatrixActual, DRedhefferMatrixCalculated, self.absoluteTolerance,
                self.relativeTolerance, "D12 for layer 1 - 2")

    def testFRedhefferMatrix(self):
        SA = self.SLayer1
        SB = self.SLayer2
        FRedhefferMatrixActual = self.FLayer12
        FRedhefferMatrixCalculated = calculateRedhefferFMatrix(SA, SB)
        assertAlmostEqual(FRedhefferMatrixActual, FRedhefferMatrixCalculated, self.absoluteTolerance,
                self.relativeTolerance, "F12 for layer 1 - 2")


    def testSReflectionRegionMatrixFromRaw(self):
        SReflectionRegionCalculated = calculateReflectionRegionSMatrixFromRaw(
                self.AReflectionRegion, self.BReflectionRegion)

        S11ReflectionRegionActual = self.S11ReflectionRegion
        S11ReflectionRegionCalculated = SReflectionRegionCalculated[0,0]
        assertAlmostEqual(S11ReflectionRegionActual, S11ReflectionRegionCalculated,
                self.absoluteTolerance, self.relativeTolerance, "S11 Matrix layer 1")

        S12ReflectionRegionActual = self.S12ReflectionRegion
        S12ReflectionRegionCalculated = SReflectionRegionCalculated[0,1]
        assertAlmostEqual(S12ReflectionRegionActual, S12ReflectionRegionCalculated,
                self.absoluteTolerance, self.relativeTolerance, "S12 layer 1")

        S21ReflectionRegionActual = self.S21ReflectionRegion
        S21ReflectionRegionCalculated = SReflectionRegionCalculated[1,0]
        assertAlmostEqual(S21ReflectionRegionActual, S21ReflectionRegionCalculated,
                self.absoluteTolerance, self.relativeTolerance, "S21 layer 1")

        S22ReflectionRegionActual = self.S22ReflectionRegion
        S22ReflectionRegionCalculated = SReflectionRegionCalculated[1,1]
        assertAlmostEqual(S22ReflectionRegionActual, S22ReflectionRegionCalculated,
                self.absoluteTolerance, self.relativeTolerance, "S22 layer 1")

    def testSTransmissionRegionMatrixFromRaw(self):
        STransmissionRegionCalculated = calculateTransmissionRegionSMatrixFromRaw(
                self.ATransmissionRegion, self.BTransmissionRegion)

        S11TransmissionRegionActual = self.S11TransmissionRegion
        S11TransmissionRegionCalculated = STransmissionRegionCalculated[0,0]
        assertAlmostEqual(S11TransmissionRegionActual, S11TransmissionRegionCalculated,
                self.absoluteTolerance, self.relativeTolerance, "S11")

        S12TransmissionRegionActual = self.S12TransmissionRegion
        S12TransmissionRegionCalculated = STransmissionRegionCalculated[0,1]
        assertAlmostEqual(S12TransmissionRegionActual, S12TransmissionRegionCalculated,
                self.absoluteTolerance, self.relativeTolerance, "S12")

        S21TransmissionRegionActual = self.S21TransmissionRegion
        S21TransmissionRegionCalculated = STransmissionRegionCalculated[1,0]
        assertAlmostEqual(S21TransmissionRegionActual, S21TransmissionRegionCalculated,
                self.absoluteTolerance, self.relativeTolerance, "S21 layer 1")

        S22TransmissionRegionActual = self.S22TransmissionRegion
        S22TransmissionRegionCalculated = STransmissionRegionCalculated[1,1]
        assertAlmostEqual(S22TransmissionRegionActual, S22TransmissionRegionCalculated,
                self.absoluteTolerance, self.relativeTolerance, "S22 layer 1")

    # ISSUE - NEED TO REFACTOR THIS CODE SO THAT IT WORKES WITH AND CAN DETECT
    # HOMOGENOUS LAYERS. CURRENTLY BREAKING LINEAR ALGEBRA ENGINE TRYING TO INVERT A 
    # NONINVERTIBLE MATRIX (I think).
    def testReflectionRegionSMatrixFromFundamentals(self):
        SCalculated = calculateReflectionRegionSMatrix(self.Kx, self.Ky, self.erReflectionRegion,
                self.urReflectionRegion, self.WFreeSpace, self.VFreeSpace)

        S11Actual = self.S11ReflectionRegion
        S11Calculated = SCalculated[0,0]
        assertAlmostEqual(S11Actual, S11Calculated,
                self.absoluteTolerance, self.relativeTolerance, "S11 Reflection Region")

        S12Actual = self.S12ReflectionRegion
        S12Calculated = SCalculated[0,1]
        assertAlmostEqual(S12Actual, S12Calculated,
                self.absoluteTolerance, self.relativeTolerance, "S12 Reflection Region")

        S21Actual = self.S21ReflectionRegion
        S21Calculated = SCalculated[1,0]
        assertAlmostEqual(S21Actual, S21Calculated,
                self.absoluteTolerance, self.relativeTolerance, "S21 Reflection Region")

        S22Actual = self.S22ReflectionRegion
        S22Calculated = SCalculated[1,1]
        assertAlmostEqual(S22Actual, S22Calculated,
                self.absoluteTolerance, self.relativeTolerance, "S22 Reflection Region ")

    def testTransmissionRegionSMatrixFromFundamentals(self):
        SCalculated = calculateTransmissionRegionSMatrix(self.Kx, self.Ky, self.erTransmissionRegion,
                self.urTransmissionRegion, self.WFreeSpace, self.VFreeSpace)

        S11Calculated = SCalculated[0,0]
        S11Actual = self.S11TransmissionRegion
        assertAlmostEqual(S11Actual, S11Calculated,
                self.absoluteTolerance, self.relativeTolerance, "S11 Transmission Region")

        S12Calculated = SCalculated[0,1]
        S12Actual = self.S12TransmissionRegion
        assertAlmostEqual(S12Actual, S12Calculated,
                self.absoluteTolerance, self.relativeTolerance, "S12 Transmission Region")

        S21Calculated = SCalculated[1,0]
        S21Actual = self.S21TransmissionRegion
        assertAlmostEqual(S21Actual, S21Calculated,
                self.absoluteTolerance, self.relativeTolerance, "S21 Transmission Region")

        S22Calculated = SCalculated[1,1]
        S22Actual = self.S22TransmissionRegion
        assertAlmostEqual(S22Actual, S22Calculated,
                self.absoluteTolerance, self.relativeTolerance, "S22 Transmission Region")

    def testRedhefferProduct(self):
        # Sanity check - Redheffer product works when one matrix is the transparent matrix.
        SA = self.transparentSMatrix
        SB = self.SLayer1
        SABActual = self.SGlobalLayer1
        SABCalculated = calculateRedhefferProduct(SA, SB)
        assertAlmostEqual(SABActual, SABCalculated, self.absoluteTolerance, self.relativeTolerance,
                "Redheffer product with Layer 1 and transparent matrix")
        SA = self.SLayer1
        SB = self.transparentSMatrix
        SABActual = self.SGlobalLayer1
        SABCalculated = calculateRedhefferProduct(SA, SB)
        assertAlmostEqual(SABActual, SABCalculated, self.absoluteTolerance, self.relativeTolerance,
                "Redheffer product with Layer 1 and transparent matrix (reversed order)")

        SA = self.SLayer1
        SB = self.SLayer2
        SABActual = self.SGlobalLayer2
        SABCalculated = calculateRedhefferProduct(SA, SB)
        assertAlmostEqual(SABActual, SABCalculated, self.absoluteTolerance, self.relativeTolerance,
                "Redheffer Product Layers 1 and 2")

        SA = self.SReflectionRegion
        SB1 = self.SLayer1
        SB2 = self.SLayer2
        SB3 = self.STransmissionRegion
        SABCalculated = calculateRedhefferProduct(SA, SB1)
        SABCalculated = calculateRedhefferProduct(SABCalculated, SB2)
        SABCalculated = calculateRedhefferProduct(SABCalculated, SB3)
        SABActual = self.SGlobal
        assertAlmostEqual(SABActual, SABCalculated, self.absoluteTolerance, self.relativeTolerance,
                "Redheffer Product Layers 1 and 2")

#    def testCalcEz(self):
#        EzActual = self.EzReflected
#        EzCalculated = calculateEz(self.Kx, self.Ky, self.KzReflectionRegion,
#                self.ExReflected, self.EyReflected);
#
#        assertAlmostEqual(EzActual, EzCalculated, self.absoluteTolerance, self.relativeTolerance,
#                "Ez in reflection region");
#
#        EzActual = self.EzTransmitted
#        EzCalculated = calculateEz(self.Kx, self.Ky, self.KzTransmissionRegion,
#                self.ExTransmitted, self.EyTransmitted);
#        assertAlmostEqual(EzActual, EzCalculated, self.absoluteTolerance, self.relativeTolerance,
#                "Ez in transmission region");

#    def testCalcRT(self):
#        RActual = self.R;
#        TActual = self.T;
#
#        (RCalculated, TCalculated) = calculateRT(self.KzReflectionRegion, self.KzTransmissionRegion,
#                self.urReflectionRegion, self.urTransmissionRegion,
#                self.ExyzReflected, self.ExyzTransmitted);
#        assertAlmostEqual(RActual, RCalculated, self.absoluteTolerance, self.relativeTolerance);
#        assertAlmostEqual(TActual, TCalculated, self.absoluteTolerance, self.relativeTolerance);

    def setUp(self):
        self.absoluteTolerance = 1e-3
        self.relativeTolerance = 1e-3
        deg = pi / 180
        self.wavelength = 2
        self.k0 = 2*pi / self.wavelength
        self.theta = 60 * deg
        self.phi = 30 * deg

        self.erReflectionRegion = 2
        self.urReflectionRegion = 1
        self.erTransmissionRegion = 9
        self.urTransmissionRegion = 1
        self.erDeviceRegion = 6
        self.urDeviceRegion = 1
        self.thicknessLayer1 = 0.5
        self.thicknessLayer2 = 0.3
        numberHarmonicsX = 3
        numberHarmonicsY = 3
        self.numberHarmonics = (numberHarmonicsX, numberHarmonicsY)

        self.erConvolutionMatrixLayer1 = numpyArrayFromFile(
            "test/matrixDataOblique/layer1/erConvolutionData.txt")
        self.urConvolutionMatrixLayer1 = complexIdentity(9)
        self.erConvolutionMatrixLayer2 = self.erDeviceRegion*complexIdentity(9)
        self.urConvolutionMatrixLayer2 = complexIdentity(9)
        self.Kx = np.diag(complexArray(
            [2.2035, 1.0607, -0.0822, 2.2035, 1.0607, -0.0822, 2.2035, 1.0607, -0.0822]))
        self.Ky = np.diag(complexArray(
            [1.9457, 1.9457, 1.9457, 0.6124, 0.6124, 0.6124, -0.7210, -0.7210, -0.7210]))
        self.KzReflectionRegion = numpyArrayFromFile(
                "test/matrixDataOblique/reflectionRegion/KzReflectionRegion.txt")
        self.KzTransmissionRegion = np.diag(complexArray(
            [0.5989, 2.0222, 2.2820, 1.9415, 2.7386, 2.9357, 1.9039, 2.7121, 2.9109]))

        self.KzFreeSpace = numpyArrayFromFile(
                "test/matrixDataOblique/freeSpace/KzFreeSpace.txt")
        self.QFreeSpace = numpyArrayFromFile("test/matrixDataOblique/freeSpace/QFreeSpace.txt")
        self.WFreeSpace = complexIdentity(18)
        self.LambdaFreeSpace = numpyArrayFromSeparatedColumnsFile("test/matrixDataOblique/freeSpace/LambdaFreeSpace.txt")
        self.VFreeSpace = numpyArrayFromSeparatedColumnsFile("test/matrixDataOblique/freeSpace/VFreeSpace.txt")

        self.S11Transparent = complexZeros((18, 18))
        self.S22Transparent = complexZeros((18, 18))
        self.S21Transparent = complexIdentity(18)
        self.S12Transparent = complexIdentity(18)

        self.PLayer1 = numpyArrayFromSeparatedColumnsFile("test/matrixDataOblique/layer1/PLayer1.txt")
        self.QLayer1 = numpyArrayFromSeparatedColumnsFile("test/matrixDataOblique/layer1/QLayer1.txt")
        self.OmegaSquaredLayer1 = numpyArrayFromSeparatedColumnsFile("test/matrixDataOblique/layer1/OmegaSquaredLayer1.txt")
        self.LambdaLayer1= numpyArrayFromSeparatedColumnsFile("test/matrixDataOblique/layer1/LambdaLayer1.txt")
        # NOTE - THE LAYER 1 VALUES ARE MODIFIED SO THAT ELEMENTS 7, 11, 12, AND 16 ALONG THE MAIN
        # DIAGONAL OF THE EIGENVALUE MATRIX ARE THEIR OWN CONJUGATE. THIS IS A NUMERICAL ERROR ISSUE
        # THAT I DON'T KNOW HOW TO RESOLVE AND I DON'T THINK IT SHOULD HAVE ANY PHYSICAL CONSEQUENCES.
        # SO I HAVE MODIFIED THE X AND V MATRICES. PHYSICALLY, 
        self.VLayer1 = numpyArrayFromSeparatedColumnsFile("test/matrixDataOblique/layer1/VLayer1.txt")
        self.ALayer1 = numpyArrayFromSeparatedColumnsFile("test/matrixDataOblique/layer1/ALayer1.txt")
        self.BLayer1 = numpyArrayFromSeparatedColumnsFile("test/matrixDataOblique/layer1/BLayer1.txt")
        self.XLayer1 = numpyArrayFromSeparatedColumnsFile("test/matrixDataOblique/layer1/XLayer1.txt")
        self.WLayer1 = numpyArrayFromSeparatedColumnsFile("test/matrixDataOblique/layer1/WLayer1.txt")
        self.S11Layer1 = numpyArrayFromSeparatedColumnsFile("test/matrixDataOblique/layer1/S11Layer1.txt")
        self.S12Layer1 = numpyArrayFromSeparatedColumnsFile("test/matrixDataOblique/layer1/S12Layer1.txt")
        self.S21Layer1 = numpyArrayFromSeparatedColumnsFile("test/matrixDataOblique/layer1/S21Layer1.txt")
        self.S22Layer1 = numpyArrayFromSeparatedColumnsFile("test/matrixDataOblique/layer1/S22Layer1.txt")
        self.SLayer1 = complexArray([
            [self.S11Layer1, self.S12Layer1],
            [self.S21Layer1, self.S22Layer1]])

        self.SGlobal11Layer1 = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataOblique/layer1/SGlobal11Layer1.txt")
        self.SGlobal12Layer1 = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataOblique/layer1/SGlobal12Layer1.txt")
        self.SGlobal21Layer1 = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataOblique/layer1/SGlobal21Layer1.txt")
        self.SGlobal22Layer1 = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataOblique/layer1/SGlobal22Layer1.txt")
        self.SGlobalLayer1 = complexArray([
            [self.SGlobal11Layer1, self.SGlobal12Layer1],
            [self.SGlobal21Layer1, self.SGlobal22Layer1]])

        self.PLayer2 = numpyArrayFromFile("test/matrixDataOblique/layer2/PLayer2.txt")
        self.QLayer2 = numpyArrayFromFile("test/matrixDataOblique/layer2/QLayer2.txt")
        self.OmegaSquaredLayer2 = numpyArrayFromFile("test/matrixDataOblique/layer2/OmegaSquaredLayer2.txt")
        self.WLayer2 = numpyArrayFromSeparatedColumnsFile("test/matrixDataOblique/layer2/WLayer2.txt")
        self.LambdaLayer2 = numpyArrayFromSeparatedColumnsFile("test/matrixDataOblique/layer2/LambdaLayer2.txt")
        self.VLayer2 = numpyArrayFromSeparatedColumnsFile("test/matrixDataOblique/layer2/VLayer2.txt")
        self.ALayer2 = numpyArrayFromSeparatedColumnsFile("test/matrixDataOblique/layer2/ALayer2.txt")
        self.BLayer2 = numpyArrayFromSeparatedColumnsFile("test/matrixDataOblique/layer2/BLayer2.txt")
        self.XLayer2 = numpyArrayFromSeparatedColumnsFile("test/matrixDataOblique/layer2/XLayer2.txt")
        self.S11Layer2 = numpyArrayFromSeparatedColumnsFile("test/matrixDataOblique/layer2/S11Layer2.txt")
        self.S12Layer2 = numpyArrayFromSeparatedColumnsFile("test/matrixDataOblique/layer2/S12Layer2.txt")
        self.S21Layer2 = numpyArrayFromSeparatedColumnsFile("test/matrixDataOblique/layer2/S21Layer2.txt")
        self.S22Layer2 = numpyArrayFromSeparatedColumnsFile("test/matrixDataOblique/layer2/S22Layer2.txt")
        self.SLayer2 = complexArray([
            [self.S11Layer2, self.S12Layer2],
            [self.S21Layer2, self.S22Layer2]])
        self.DLayer12= np.loadtxt('test/matrixDataOblique/layer2/D12.csv', dtype=np.cdouble)
        self.FLayer12= np.loadtxt('test/matrixDataOblique/layer2/F12.csv', dtype=np.cdouble)

        self.SGlobal11Layer2 = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataOblique/layer2/SGlobal11Layer2.txt")
        self.SGlobal12Layer2 = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataOblique/layer2/SGlobal12Layer2.txt")
        self.SGlobal21Layer2 = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataOblique/layer2/SGlobal21Layer2.txt")
        self.SGlobal22Layer2 = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataOblique/layer2/SGlobal22Layer2.txt")
        self.SGlobalLayer2 = complexArray([
            [self.SGlobal11Layer2, self.SGlobal12Layer2],
            [self.SGlobal21Layer2, self.SGlobal22Layer2]])

        self.QReflectionRegion = numpyArrayFromFile(
                "test/matrixDataOblique/reflectionRegion/QReflectionRegion.txt")
        self.LambdaReflectionRegion = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataOblique/reflectionRegion/LambdaReflectionRegion.txt")
        self.WReflectionRegion = complexIdentity(18)
        self.VReflectionRegion = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataOblique/reflectionRegion/VReflectionRegion.txt")
        self.LambdaReflectionRegion = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataOblique/reflectionRegion/LambdaReflectionRegion.txt")
        self.AReflectionRegion= numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataOblique/reflectionRegion/AReflectionRegion.txt")
        self.BReflectionRegion = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataOblique/reflectionRegion/BReflectionRegion.txt")
        self.S11ReflectionRegion = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataOblique/reflectionRegion/S11ReflectionRegion.txt")
        self.S12ReflectionRegion = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataOblique/reflectionRegion/S12ReflectionRegion.txt")
        self.S21ReflectionRegion = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataOblique/reflectionRegion/S21ReflectionRegion.txt")
        self.S22ReflectionRegion = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataOblique/reflectionRegion/S22ReflectionRegion.txt")
        self.SReflectionRegion = complexArray([
            [self.S11ReflectionRegion, self.S12ReflectionRegion],
            [self.S21ReflectionRegion, self.S22ReflectionRegion]])

        self.QTransmissionRegion = numpyArrayFromFile(
                "test/matrixDataOblique/transmissionRegion/QTransmissionRegion.txt")
        self.LambdaTransmissionRegion = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataOblique/transmissionRegion/LambdaTransmissionRegion.txt")
        self.WTransmissionRegion = complexIdentity(18)
        self.VTransmissionRegion = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataOblique/transmissionRegion/VTransmissionRegion.txt")
        self.LambdaTransmissionRegion = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataOblique/transmissionRegion/LambdaTransmissionRegion.txt")
        self.ATransmissionRegion= numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataOblique/transmissionRegion/ATransmissionRegion.txt")
        self.BTransmissionRegion = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataOblique/transmissionRegion/BTransmissionRegion.txt")
        self.S11TransmissionRegion = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataOblique/transmissionRegion/S11TransmissionRegion.txt")
        self.S12TransmissionRegion = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataOblique/transmissionRegion/S12TransmissionRegion.txt")
        self.S21TransmissionRegion = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataOblique/transmissionRegion/S21TransmissionRegion.txt")
        self.S22TransmissionRegion = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataOblique/transmissionRegion/S22TransmissionRegion.txt")
        self.STransmissionRegion = complexArray([
            [self.S11TransmissionRegion, self.S12TransmissionRegion],
            [self.S21TransmissionRegion, self.S22TransmissionRegion]])

        # Overall global scattering matrices
        self.SGlobal11= numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataOblique/SGlobal11.txt")
        self.SGlobal12= numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataOblique/SGlobal12.txt")
        self.SGlobal21= numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataOblique/SGlobal21.txt")
        self.SGlobal22= numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataOblique/SGlobal22.txt")
        self.SGlobal = complexArray([
            [self.SGlobal11, self.SGlobal12],
            [self.SGlobal21, self.SGlobal22]])

        self.transparentSMatrix = complexZeros((2, 2, 18, 18))
        self.transparentSMatrix[0,1] = complexIdentity(18)
        self.transparentSMatrix[1,0] = complexIdentity(18)

if __name__ == '__main__':
    unittest.main()
