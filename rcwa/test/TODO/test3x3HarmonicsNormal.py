import sys
sys.path.append('core');
sys.path.append('test')

import unittest
from shorthandTest import *
from matrices import *
from fresnel import *
from convolution import generateConvolutionMatrix
from matrixParser import *


class Test3x3HarmonicsNormal(unittest.TestCase):
    def testCalculateKz(self):
        KzCalculated = calculateKz(self.Kx, self.Ky, self.erReflectionRegion, self.urReflectionRegion)
        KzActual = self.KzReflectionRegion
        assertAlmostEqual(KzActual, KzCalculated,
                self.absoluteTolerance, self.relativeTolerance, "Kz in Reflection region not correct")

        KzCalculated = calculateKz(self.Kx, self.Ky, self.erTransmissionRegion, self.urTransmissionRegion)
        KzActual = self.KzTransmissionRegion
        assertAlmostEqual(KzActual, KzCalculated,
                self.absoluteTolerance, self.relativeTolerance, "Kz in transmission region not correct")

    def testaTEM(self):
        aTEActual = complexArray([-0.39073, 0.920505])
        aTMActual = complexArray([0.50137, 0.21282])
        (aTECalculated, aTMCalculated) = aTEMGen(self.Kx, self.Ky, self.KzReflectionRegion);
        assertAlmostEqual(aTEActual, aTECalculated,
                self.absoluteTolerance, self.relativeTolerance, "Oblique incidence TE vector wrong");
        assertAlmostEqual(aTMActual, aTMCalculated,
                self.absoluteTolerance, self.relativeTolerance, "Oblique incidence TM vector wrong");

        # Important corner case where the cross product fails
        kx = 0;
        ky = 0.0001;
        kz = 1.5;
        aTEActual = complexArray([0,1,0])[0:2];
        aTMActual = complexArray([1,0,0])[0:2];
        (aTECalculated, aTMCalculated) = aTEMGen(kx, ky, kz);
        assertAlmostEqual(aTEActual, aTECalculated,
                self.absoluteTolerance, self.relativeTolerance, "Near-normal Incidence TE vector wrong");
        assertAlmostEqual(aTMActual, aTMCalculated,
                self.absoluteTolerance, self.relativeTolerance, "Near-normal incidence TM vector wrong");

    def testTransparentSMatrix(self):

        SActual = self.transparentSMatrix
        SCalculated = generateTransparentSMatrix();

        assertAlmostEqual(SActual, SCalculated,self.absoluteTolerance,self.relativeTolerance);

    def testCalculateKVector(self):
        kVectorActual = complexArray([self.Kx, self.Ky, self.KzReflectionRegion])
        kVectorCalculated = calculateKVector(self.theta, self.phi,
                self.erReflectionRegion, self.urReflectionRegion)
        assertAlmostEqual(kVectorActual, kVectorCalculated,self.absoluteTolerance,self.relativeTolerance);

    def testCalcEz(self):
        EzActual = self.EzReflected
        EzCalculated = calculateEz(self.Kx, self.Ky, self.KzReflectionRegion,
                self.ExReflected, self.EyReflected);

        assertAlmostEqual(EzActual, EzCalculated, self.absoluteTolerance, self.relativeTolerance,
                "Ez in reflection region");

        EzActual = self.EzTransmitted
        EzCalculated = calculateEz(self.Kx, self.Ky, self.KzTransmissionRegion,
                self.ExTransmitted, self.EyTransmitted);
        assertAlmostEqual(EzActual, EzCalculated, self.absoluteTolerance, self.relativeTolerance,
                "Ez in transmission region");

    def testCalcRT(self):
        RActual = self.R;
        TActual = self.T;

        (RCalculated, TCalculated) = calculateRT(self.KzReflectionRegion, self.KzTransmissionRegion,
                self.urReflectionRegion, self.urTransmissionRegion,
                self.ExyzReflected, self.ExyzTransmitted);
        assertAlmostEqual(RActual, RCalculated, self.absoluteTolerance, self.relativeTolerance);
        assertAlmostEqual(TActual, TCalculated, self.absoluteTolerance, self.relativeTolerance);
    def testPMatrix(self):
        PActual = complexArray([
            [0.212504, 0.499373],
            [-0.909798, -0.212504]])
        PCalculated = calculatePMatrix(self.Kx, self.Ky, self.erLayer1, self.urLayer1)
        assertAlmostEqual(PActual, PCalculated, self.absoluteTolerance, self.relativeTolerance,
                "P matrix layer 1");

    def testQMatrix(self):
        QActual = complexArray([[0.4250, 0.9987],[-1.8196, -0.4250]]);
        QCalculated = calculateQMatrix(self.Kx, self.Ky, self.erLayer1, self.urLayer1);
        assertAlmostEqual(QActual, QCalculated, self.absoluteTolerance, self.relativeTolerance,
                "Q matrix Layer 1");

        QActual = complexArray([[0.1417, 0.6662],[-0.9399, -0.1417]]);
        QCalculated = calculateQMatrix(self.Kx, self.Ky, self.erLayer2, self.urLayer2);
        assertAlmostEqual(QActual, QCalculated, self.absoluteTolerance, self.relativeTolerance,
                "Q matrix layer 2")

    def testOmegaMatrix(self):
        OActual = complexArray([[0 + 0.9046j, 0+0j],[0+0j,0+0.9046j]]);
        OCalculated = calculateOmegaMatrix(self.KzLayer1);
        assertAlmostEqual(OActual, OCalculated, self.absoluteTolerance, self.relativeTolerance);

        OActual = complexArray([[0 + 1.3485j, 0+0j],[0+0j,0+1.3485j]]);
        OCalculated = calculateOmegaMatrix(self.KzLayer2);
        assertAlmostEqual(OActual, OCalculated, self.absoluteTolerance, self.relativeTolerance);

    def testVMatrix(self):
        (VCalculated, W) = calculateVWXMatrices(self.Kx, self.Ky, self.KzGap,
                self.erGap, self.urGap)
        VActual = complexArray([[0 - 0.4250j, 0 - 1.1804j], [0 + 2.0013j, 0 + 0.4250j]]);
        assertAlmostEqual(VActual, VCalculated, self.absoluteTolerance, self.relativeTolerance);

        (VCalculated, W) = calculateVWXMatrices(self.Kx, self.Ky, self.KzLayer1,
                self.erLayer1, self.urLayer1)
        VActual = complexArray([[0-0.4698j,0-1.1040j],[0+2.0114j,0+0.4698j]]);
        assertAlmostEqual(VActual, VCalculated, self.absoluteTolerance, self.relativeTolerance);

        (VCalculated, W) = calculateVWXMatrices(self.Kx, self.Ky, self.KzLayer2,
                self.erLayer2, self.urLayer2)
        VActual = complexArray([[0-0.1051j,0-0.4941j],[0+0.6970j,0+0.1051j]]);
        assertAlmostEqual(VActual, VCalculated, self.absoluteTolerance, self.relativeTolerance);

        (VCalculated, W_ref) = calculateVWXMatrices(self.Kx, self.Ky, self.KzReflectionRegion,
                self.erReflectionRegion,
                self.urReflectionRegion)
        VActual = complexArray([
            [0 - 0.5017j, 0 - 0.8012j],
            [0 + 1.7702j, 0 + 0.5017j]]);
        assertAlmostEqual(VActual, VCalculated, self.absoluteTolerance, self.relativeTolerance);

    def testXMatrix(self):
        (V, W, XCalculated) = calculateVWXMatrices(self.Kx, self.Ky, self.KzLayer1,
                self.erLayer1, self.urLayer1, self.k0, self.thicknessLayer1)
        XActual = complexArray([[0.1493+0.9888j, 0+0j],[0+0j,0.1493+0.9888j]]);
        assertAlmostEqual(XActual, XCalculated, self.absoluteTolerance, self.relativeTolerance);

        (V, W, XCalculated) = calculateVWXMatrices(self.Kx, self.Ky, self.KzLayer2,
                self.erLayer2, self.urLayer2, self.k0, self.thicknessLayer2)
        XActual = complexArray([[-0.4583 - 0.8888j, 0+0j],[0+0j, -0.4583 - 0.8888j]]);
        assertAlmostEqual(XActual, XCalculated, self.absoluteTolerance, self.relativeTolerance);


    def testAMatrix(self):
        W1 = np.identity(2);
        Wg = np.identity(2);
        V1 = complexArray([[0 - 0.4698j, 0 - 1.1040j],[0 + 2.0114j, 0 + 0.4698j]]);
        Vg = complexArray([[0 - 0.4250j, 0 - 1.1804j], [0 + 2.0013j, 0 + 0.4250j]]);

        ACalculated = calculateScatteringAMatrix(W1, Wg, V1, Vg);
        AActual = self.ALayer1
        assertAlmostEqual(AActual, ACalculated, self.absoluteTolerance, self.relativeTolerance);

        W2 = complexIdentity(2);
        Wg = complexIdentity(2);
        V2 = complexArray([[0 - 0.1051j, 0 - 0.4941j],[0 + 0.6970j, 0 + 0.1051j]]);
        Vg = complexArray([[0 - 0.4250j, 0 - 1.1804j],[0 + 2.0013j, 0 + 0.4250j]]);

        ACalculated = calculateScatteringAMatrix(W2, Wg, V2, Vg);
        AActual = self.ALayer2
        assertAlmostEqual(AActual, ACalculated, self.absoluteTolerance, self.relativeTolerance);

    def testBMatrix(self):
        W1 = complexIdentity(2);
        Wg = complexIdentity(2);
        V1 = complexArray([[0 - 0.4698j, 0 - 1.1040j],[0 + 2.0114j, 0 + 0.4698j]]);
        Vg = complexArray([[0 - 0.4250j, 0 - 1.1804j], [0 + 2.0013j, 0 + 0.4250j]]);
        BCalculated = calculateScatteringBMatrix(W1, Wg, V1, Vg);
        BActual = complexArray([[-0.0049, 0.0427],[0.0427, -0.0873]]);
        assertAlmostEqual(BActual, BCalculated, self.absoluteTolerance, self.relativeTolerance);

        W2 = complexIdentity(2);
        Wg = complexIdentity(2);
        V2 = complexArray([[0 - 0.1051j, 0 - 0.4941j],[0 + 0.6970j, 0 + 0.1051j]]);
        Vg = complexArray([[0 - 0.4250j, 0 - 1.1804j],[0 + 2.0013j, 0 + 0.4250j]]);

        BCalculated = calculateScatteringBMatrix(W2, Wg, V2, Vg);
        BActual = complexArray([[-1.8324, -0.2579],[-0.2579, -1.3342]]);
        assertAlmostEqual(BActual, BCalculated, self.absoluteTolerance, self.relativeTolerance);

    def testDiMatrix(self):
        absoluteTolerance = 0.003;# D matrix has some very small values after multiple matrix mult.
        relativeTolerance = 0.1; # relative error is huge on 5e-4 values. Overriding.

        A = complexArray([[2.0049, -0.0427], [-0.0427, 2.0873]]);
        B = complexArray([[-0.0049, 0.0427], [0.0427, -0.0873]]);
        X = complexArray([[0.1493 + 0.9888j, 0+0j],[0+0j, 0.4193 + 0.9888j]]);
        DCalculated = calculateScatteringDMatrix(A, B, X);
        DActual = complexArray([[2.0057 - 0.0003j, -0.0445 + 0.0006j],[-0.0445 + 0.0006j, 2.0916 - 0.0013j]])
        assertAlmostEqual(DActual, DCalculated, absoluteTolerance, relativeTolerance);

        # LAYER 2 DATA
        # Since now we have the d-matrix to higher precision we can test it more strongly.
        A = complexArray([[3.8324, 0.2579],[0.2579, 3.3342]]);
        B = complexArray([[-1.8324, -0.2579], [-0.2579, -1.3342]]);
        X = complexArray([[-0.4583 - 0.8888j, 0+0j],[0+0j, -0.4583 - 0.8888j]]);

        DCalculated = calculateScatteringDMatrix(A, B, X);
        DActual = complexArray([[4.3436 - 0.7182j, 0.3604 - 0.1440j], [0.3604 - 0.1440j, 3.6475 - 0.4401j]]);
        assertAlmostEqual(DActual, DCalculated, self.absoluteTolerance, self.relativeTolerance);

    def testScatteringMatrixFromRaw(self):
        SMatrixLayer1Calculated = calculateInternalSMatrixFromRaw(self.ALayer1, self.BLayer1,
                self.XLayer1, self.DLayer1)
        SMatrixLayer2Calculated = calculateInternalSMatrixFromRaw(self.ALayer2, self.BLayer2,
                self.XLayer2, self.DLayer2)

        S11Actual = self.S11Layer1
        S11Calculated = SMatrixLayer1Calculated[0,0];
        assertAlmostEqual(S11Actual, S11Calculated, self.absoluteTolerance, self.relativeTolerance,
                "S11 for Layer 1");
        S11Calculated = SMatrixLayer2Calculated[0,0];
        S11Actual = self.S11Layer2
        assertAlmostEqual(S11Actual, S11Calculated, self.absoluteTolerance, self.relativeTolerance,
                "S11 for Layer 2");

        S12Actual = self.S12Layer1
        S12Calculated = SMatrixLayer1Calculated[0,1];
        assertAlmostEqual(S12Actual, S12Calculated, self.absoluteTolerance, self.relativeTolerance,
                "S12 for Layer 1");
        S12Actual = self.S12Layer2
        S12Calculated = SMatrixLayer2Calculated[0,1];
        assertAlmostEqual(S12Actual, S12Calculated, self.absoluteTolerance, self.relativeTolerance,
                "S12 for Layer 2");

        S21Actual = self.S21Layer1
        S21Calculated = SMatrixLayer1Calculated[1,0];
        assertAlmostEqual(S21Actual, S21Calculated, self.absoluteTolerance, self.relativeTolerance,
                "S21 for Layer 1");
        S21Actual = self.S21Layer2
        S21Calculated = SMatrixLayer2Calculated[1,0];
        assertAlmostEqual(S21Actual, S21Calculated, self.absoluteTolerance, self.relativeTolerance,
                "S21 for Layer 2");

        S22Actual = self.S22Layer1
        S22Calculated = SMatrixLayer1Calculated[1,1];
        assertAlmostEqual(S22Actual, S22Calculated, self.absoluteTolerance, self.relativeTolerance,
                "S22 for Layer 1");
        S22Actual = self.S22Layer2
        S22Calculated = SMatrixLayer2Calculated[1,1];
        assertAlmostEqual(S22Actual, S22Calculated, self.absoluteTolerance, self.relativeTolerance,
                "S22 for Layer 2");

    def testDRedhefferMatrix(self):
        SA = self.transparentSMatrix
        SB = self.SLayer1
        DRedhefferMatrixActual = complexArray([[1,0],[0,1]])
        DRedhefferMatrixCalculated = calculateRedhefferDMatrix(SA, SB)
        assertAlmostEqual(DRedhefferMatrixActual, DRedhefferMatrixCalculated, self.absoluteTolerance,
                self.relativeTolerance, "Layer 1 D matrix")

        SA = self.SLayer1
        SB = self.SLayer2
        DRedhefferMatrixActual = complexArray([
            [0.1506 + 0.9886j, -0.0163 - 0.0190j],
            [-0.0163 - 0.0190j, 0.1822 + 1.0253j]]);
        DRedhefferMatrixCalculated = calculateRedhefferDMatrix(SA, SB)
        assertAlmostEqual(DRedhefferMatrixActual, DRedhefferMatrixCalculated, self.absoluteTolerance,
                self.relativeTolerance, "Layer 2 D matrix")

    def testFRedhefferMatrix(self):
        SA = self.transparentSMatrix
        SB = self.SLayer1
        FRedhefferMatrixActual = complexArray([
            [0.1490 + 0.9880j, 0.0005 + 0.0017j],
            [0.0005 + 0.0017j, 0.148 + 0.9848j]]);
        FRedhefferMatrixCalculated = calculateRedhefferFMatrix(SA, SB)
        assertAlmostEqual(FRedhefferMatrixActual, FRedhefferMatrixCalculated, self.absoluteTolerance,
                self.relativeTolerance, "Layer 1 F matrix")

        SA = self.SLayer1
        SB = self.SLayer2
        FRedhefferMatrixActual = complexArray([
            [-0.2117 - 0.6413j, 0.0471 + 0.0518j],
            [0.0471 + 0.0518j, -0.3027 - 0.7414j]]);
        FRedhefferMatrixCalculated = calculateRedhefferFMatrix(SA, SB)
        assertAlmostEqual(FRedhefferMatrixActual, FRedhefferMatrixCalculated, self.absoluteTolerance,
                self.relativeTolerance, "Layer 2 F matrix")

    def testRedhefferProduct(self):
        SA = self.transparentSMatrix
        SB = self.SLayer1
        SABActual = self.SLayer1
        SABCalculated = calculateRedhefferProduct(SA, SB)
        assertAlmostEqual(SABActual, SABCalculated, self.absoluteTolerance, self.relativeTolerance,
                "Redheffer product with Layer 1 and transparent matrix")

        SA = self.SLayer1
        SB = self.SLayer2
        SABActual = complexZeros((2,2,2,2));
        SABActual[0,0] = complexArray([
            [-0.5961 + 0.4214j, -0.0840 + 0.0085j],
            [-0.0840 + 0.0085j, -0.4339 + 0.4051j]]);
        SABActual[0,1] = complexArray([
            [0.6020 - 0.3046j, -0.0431 + 0.0534j],
            [-0.0431 + 0.0534j, 0.6852 - 0.4078j]]);
        SABActual[1,0] = complexArray([
            [0.6020 - 0.3046j, -0.0431 + 0.0534j],
            [-0.0431 + 0.0534j, 0.6852 - 0.4078j]]);
        SABActual[1,1] = complexArray([
            [0.6971 - 0.2216j, 0.0672 - 0.0211j],
            [0.0672 - 0.0211j, 0.5673 - 0.1808j]]);
        SABCalculated = calculateRedhefferProduct(SA, SB)
        assertAlmostEqual(SABActual, SABCalculated, self.absoluteTolerance, self.relativeTolerance,
                "Redheffer product with Layer 1 Layer 2")

    def testSMatrixFromFundamentals(self):
        SiActual = self.SLayer1
        SiCalculated = calculateInternalSMatrix(self.Kx, self.Ky, self.erLayer1, self.urLayer1,
                self.k0, self.thicknessLayer1, self.WGap, self.VGap)
        assertAlmostEqual(SiActual, SiCalculated,
                self.absoluteTolerance, self.relativeTolerance, "S Matrix layer 1")

        SiActual = self.SLayer2
        SiCalculated = calculateInternalSMatrix(self.Kx, self.Ky, self.erLayer2, self.urLayer2,
                self.k0, self.thicknessLayer2, self.WGap, self.VGap)

    def testSReflectionRegionMatrixFromRaw(self):
        absoluteTolerance = 0.007
        relativeTolerance = 0.03
        AReflectionRegion = complexArray([
            [1.86002, 0.113614],
            [0.115376, 1.64547]]);
        BReflectionRegion = complexArray([
            [0.139976, -0.113614],
            [-0.115376, 0.354529]]);
        SReflectionRegionActual = self.SReflectionRegion
        SReflectionRegionCalculated = calculateReflectionRegionSMatrixFromRaw(
                AReflectionRegion, BReflectionRegion)
        assertAlmostEqual(SReflectionRegionActual, SReflectionRegionCalculated,
                absoluteTolerance, relativeTolerance, "S Matrix layer 1")

    def testSTransmissionRegionMatrixFromRaw(self):
        absoluteTolerance = 0.007
        relativeTolerance = 0.03
        ATransmissionRegion = complexArray([
            [1.660774, -0.0652574],
            [-0.06525902, 1.786816]]);
        BTransmissionRegion = complexArray([
            [0.339226, 0.0652574],
            [0.06525902, 0.21318382]]);
        STransmissionRegionActual = self.STransmissionRegion
        STransmissionRegionCalculated = calculateTransmissionRegionSMatrixFromRaw(
                ATransmissionRegion, BTransmissionRegion)
        assertAlmostEqual(STransmissionRegionActual, STransmissionRegionCalculated,
                absoluteTolerance, relativeTolerance, "S Matrix layer 1")

    def testReflectionRegionSMatrixFromFundamentals(self):
        SActual = self.SReflectionRegion
        SCalculated = calculateReflectionRegionSMatrix(self.Kx, self.Ky, self.erReflectionRegion,
                self.urReflectionRegion, self.WGap, self.VGap)
        assertAlmostEqual(SActual, SCalculated,
                self.absoluteTolerance, self.relativeTolerance, "S Matrix layer 1")

    def testTransmissionRegionSMatrixFromFundamentals(self):
        SActual = self.STransmissionRegion
        SCalculated = calculateTransmissionRegionSMatrix(self.Kx, self.Ky, self.erTransmissionRegion,
                self.urTransmissionRegion, self.WGap, self.VGap)
        assertAlmostEqual(SActual, SCalculated,
                self.absoluteTolerance, self.relativeTolerance, "S Matrix layer 1")


    def setUp(self):
        self.setupData1x1Harmonics() # data is equal to the old stuff, with some changes.
        numberHarmonicsX = 3
        numberHarmonicsY = 3
        self.numberHarmonics = (numberHarmonicsX, numberHarmonicsY)

        self.erConvolutionMatrixLayer1 = numpyArrayFromFile(
            "test/matrixDataNormal/layer1/erConvolutionData.txt")
        self.urConvolutionMatrixLayer1 = complexIdentity(9)
        self.erConvolutionMatrixLayer2 = 6*complexIdentity(9)
        self.urConvolutionMatrixLayer2 = complexIdentity(9)
        self.Kx = np.diag(1.1429*complexArray([1,0,-1,1,0,-1,1,0,-1]))
        self.Ky = np.diag(1.3333*complexArray([1,1,1,0,0,0,-1,-1,-1]))
        self.KzReflectionRegion = numpyArrayFromFile(
                "test/matrixDataNormal/reflectionRegion/KzReflectionRegion.txt")
        self.KzTransmissionRegion = np.diag(complexArray([2.4324, 2.6874, 2.4324, 2.7738, 3.0, 2.7738,
            2.4324, 2.6874, 2.4324]))

        self.KzFreeSpace = numpyArrayFromFile("test/matrixDataNormal/freeSpace/KzFreeSpace.txt")
        self.QFreeSpace = numpyArrayFromFile("test/matrixDataNormal/freeSpace/QFreeSpace.txt")
        self.WFreeSpace = complexIdentity(18)
        self.LamFreeSpace = numpyArrayFromFile("test/matrixDataNormal/freeSpace/LamFreeSpace.txt")
        self.VFreeSpace = numpyArrayFromFile("test/matrixDataNormal/freeSpace/VFreeSpace.txt")

        self.S11Transparent = complexZeros((18, 18))
        self.S22Transparent = complexZeros((18, 18))
        self.S21Transparent = complexIdentity(18)
        self.S12Transparent = complexIdentity(18)

        self.PLayer1 = numpyArrayFromFile("test/matrixDataNormal/layer1/PLayer1.txt")
        self.QLayer1 = numpyArrayFromFile("test/matrixDataNormal/layer1/QLayer1.txt")
        self.OmegaSquaredLayer1 = numpyArrayFromFile("test/matrixDataNormal/layer1/OmegaSquaredLayer1.txt")
        self.LambdaLayer1= numpyArrayFromFile("test/matrixDataNormal/layer1/LambdaLayer1.txt")
        self.VLayer1 = numpyArrayFromFile("test/matrixDataNormal/layer1/VLayer1.txt")
        self.ALayer1 = numpyArrayFromFile("test/matrixDataNormal/layer1/ALayer1.txt")
        self.BLayer1 = numpyArrayFromFile("test/matrixDataNormal/layer1/BLayer1.txt")
        self.XLayer1 = numpyArrayFromFile("test/matrixDataNormal/layer1/XLayer1.txt")
        self.S11Layer1 = numpyArrayFromFile("test/matrixDataNormal/layer1/S11Layer1.txt")
        self.S12Layer1 = numpyArrayFromFile("test/matrixDataNormal/layer1/S12Layer1.txt")
        self.S21Layer1 = numpyArrayFromSeparatedColumnsFile("test/matrixDataNormal/layer1/S21Layer1.txt")
        self.S22Layer1 = numpyArrayFromSeparatedColumnsFile("test/matrixDataNormal/layer1/S22Layer1.txt")
        self.SGlobal11Layer1 = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataNormal/layer1/SGlobal11Layer1.txt")
        self.SGlobal12Layer1 = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataNormal/layer1/SGlobal12Layer1.txt")
        self.SGlobal21Layer1 = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataNormal/layer1/SGlobal21Layer1.txt")
        self.SGlobal22Layer1 = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataNormal/layer1/SGlobal22Layer1.txt")

        self.PLayer2 = numpyArrayFromFile("test/matrixDataNormal/layer2/PLayer2.txt")
        self.QLayer2 = numpyArrayFromFile("test/matrixDataNormal/layer2/QLayer2.txt")
        self.OmegaSquaredLayer2 = numpyArrayFromFile("test/matrixDataNormal/layer2/OmegaSquaredLayer2.txt")
        self.WLayer2 = numpyArrayFromSeparatedColumnsFile("test/matrixDataNormal/layer2/WLayer2.txt")
        self.LambdaLayer2 = numpyArrayFromSeparatedColumnsFile("test/matrixDataNormal/layer2/LambdaLayer2.txt")
        self.VLayer2 = numpyArrayFromSeparatedColumnsFile("test/matrixDataNormal/layer2/VLayer2.txt")
        self.ALayer2 = numpyArrayFromSeparatedColumnsFile("test/matrixDataNormal/layer2/ALayer2.txt")
        self.BLayer2 = numpyArrayFromSeparatedColumnsFile("test/matrixDataNormal/layer2/BLayer2.txt")
        self.XLayer2 = numpyArrayFromSeparatedColumnsFile("test/matrixDataNormal/layer2/XLayer2.txt")
        self.S11Layer2 = numpyArrayFromSeparatedColumnsFile("test/matrixDataNormal/layer2/S11Layer2.txt")
        self.S12Layer2 = numpyArrayFromSeparatedColumnsFile("test/matrixDataNormal/layer2/S12Layer2.txt")
        self.S21Layer2 = numpyArrayFromSeparatedColumnsFile("test/matrixDataNormal/layer2/S21Layer2.txt")
        self.S22Layer2 = numpyArrayFromSeparatedColumnsFile("test/matrixDataNormal/layer2/S22Layer2.txt")
        self.SGlobal11Layer2 = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataNormal/layer2/SGlobal11Layer2.txt")
        self.SGlobal12Layer2 = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataNormal/layer2/SGlobal12Layer2.txt")
        self.SGlobal21Layer2 = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataNormal/layer2/SGlobal21Layer2.txt")
        self.SGlobal22Layer2 = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataNormal/layer2/SGlobal22Layer2.txt")

        self.QReflectionRegion = numpyArrayFromFile(
                "test/matrixDataNormal/reflectionRegion/QReflectionRegion.txt")
        self.LambdaReflectionRegion = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataNormal/reflectionRegion/LambdaReflectionRegion.txt")
        self.WReflectionRegion = complexIdentity(18)
        self.VReflectionRegion = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataNormal/reflectionRegion/VReflectionRegion.txt")
        self.LambdaReflectionRegion = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataNormal/reflectionRegion/LambdaReflectionRegion.txt")
        self.AReflectionRegion= numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataNormal/reflectionRegion/AReflectionRegion.txt")
        self.BReflectionRegion = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataNormal/reflectionRegion/BReflectionRegion.txt")
        self.S11ReflectionRegion = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataNormal/reflectionRegion/S11ReflectionRegion.txt")
        self.S12ReflectionRegion = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataNormal/reflectionRegion/S12ReflectionRegion.txt")
        self.S21ReflectionRegion = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataNormal/reflectionRegion/S21ReflectionRegion.txt")
        self.S22ReflectionRegion = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataNormal/reflectionRegion/S22ReflectionRegion.txt")

        self.QTransmissionRegion = numpyArrayFromFile(
                "test/matrixDataNormal/transmissionRegion/QTransmissionRegion.txt")
        self.LambdaTransmissionRegion = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataNormal/transmissionRegion/LambdaTransmissionRegion.txt")
        self.WTransmissionRegion = complexIdentity(18)
        self.VTransmissionRegion = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataNormal/transmissionRegion/VTransmissionRegion.txt")
        self.LambdaTransmissionRegion = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataNormal/transmissionRegion/LambdaTransmissionRegion.txt")
        self.ATransmissionRegion= numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataNormal/transmissionRegion/ATransmissionRegion.txt")
        self.BTransmissionRegion = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataNormal/transmissionRegion/BTransmissionRegion.txt")
        self.S11TransmissionRegion = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataNormal/transmissionRegion/S11TransmissionRegion.txt")
        self.S12TransmissionRegion = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataNormal/transmissionRegion/S12TransmissionRegion.txt")
        self.S21TransmissionRegion = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataNormal/transmissionRegion/S21TransmissionRegion.txt")
        self.S22TransmissionRegion = numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataNormal/transmissionRegion/S22TransmissionRegion.txt")

        self.SGlobal11= numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataNormal/SGlobal11.txt")
        self.SGlobal12= numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataNormal/SGlobal12.txt")
        self.SGlobal21= numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataNormal/SGlobal21.txt")
        self.SGlobal22= numpyArrayFromSeparatedColumnsFile(
                "test/matrixDataNormal/SGlobal22.txt")
        self.delta = complexArray([0,0,0,0,1,0,0,0,0])
        # TODO: add final csrc, esrc, cref, reflection coefficients, and so on
if __name__ == '__main__':
    unittest.main()
