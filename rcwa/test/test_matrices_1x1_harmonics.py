# Tests the matrices.py file, which is responsible for the creation and manipulation of scattering matrices
import unittest
from rcwa.testing import *
from rcwa.utils import k_vector
from rcwa.shorthand import *
from rcwa import Source, Layer, LayerStack
import numpy as np
from rcwa.matrices import *

class Test1x1Harmonic(unittest.TestCase):
    def testCalculateKz(self):
        layer = self.layerStack.incident_layer
        KzCalculated = layer.Kz_forward()
        KzActual = np.float64(self.KzReflectionRegion)
        assert_almost_equal(KzActual, KzCalculated,
                            self.absoluteTolerance, self.relativeTolerance, "Kz in Reflection region not correct")

        layer = self.layerStack.transmission_layer
        KzCalculated = layer.Kz_forward()
        KzActual = self.KzTransmissionRegion
        assert_almost_equal(KzActual, KzCalculated,
                            self.absoluteTolerance, self.relativeTolerance, "Kz in transmission region not correct")

    def testaTEM(self):
        aTEActual = complexArray([-0.39073, 0.920505, 0])
        aTMActual = complexArray([0.50137, 0.21282, -0.838671])
        (aTECalculated, aTMCalculated) = self.source.aTE, self.source.aTM
        assert_almost_equal(aTEActual, aTECalculated,
                            self.absoluteTolerance, self.relativeTolerance, "Oblique incidence TE vector wrong");
        assert_almost_equal(aTMActual, aTMCalculated,
                            self.absoluteTolerance, self.relativeTolerance, "Oblique incidence TM vector wrong");

        ATEMActual = np.vstack((aTEActual, aTMActual))
        ATEMCalculated = self.source.ATEM
        assert_almost_equal(ATEMActual, ATEMCalculated,
                            self.absoluteTolerance, self.relativeTolerance, "Oblique incidence TM vector wrong");


        # Important corner case where the cross product fails
        aTEActual = complexArray([0,1,0]);
        aTMActual = complexArray([1,0,0]);
        zeroAngleSource = Source()
        (aTECalculated, aTMCalculated) = zeroAngleSource.aTE, zeroAngleSource.aTM
        assert_almost_equal(aTEActual, aTECalculated,
                            self.absoluteTolerance, self.relativeTolerance, "Near-normal Incidence TE vector wrong");
        assert_almost_equal(aTMActual, aTMCalculated,
                            self.absoluteTolerance, self.relativeTolerance, "Near-normal incidence TM vector wrong");

    def testTransparentSMatrix(self):

        SActual = self.transparentSMatrix
        SCalculated = S_matrix_transparent((2, 2));

        assert_almost_equal(SActual, SCalculated, self.absoluteTolerance, self.relativeTolerance);

    def testCalculateKVector(self):
        kVectorActual = complexArray([self.Kx, self.Ky, self.KzReflectionRegion])
        kVectorCalculated = k_vector(self.source, self.layerStack.incident_layer, normalize=True)
        assert_almost_equal(kVectorActual, kVectorCalculated, self.absoluteTolerance, self.relativeTolerance);

    def testCalcEz(self):
        EzActual = self.EzReflected
        EzCalculated = calculateEz(self.Kx, self.Ky, self.KzReflectionRegion,
                self.ExReflected, self.EyReflected);

        assert_almost_equal(EzActual, EzCalculated, self.absoluteTolerance, self.relativeTolerance,
                "Ez in reflection region");

        EzActual = self.EzTransmitted
        EzCalculated = calculateEz(self.Kx, self.Ky, self.KzTransmissionRegion,
                self.ExTransmitted, self.EyTransmitted);
        assert_almost_equal(EzActual, EzCalculated, self.absoluteTolerance, self.relativeTolerance,
                "Ez in transmission region");

    def testCalcRT(self):
        RActual = self.R;
        TActual = self.T;

        (RCalculated, TCalculated) = calculateRT(self.KzReflectionRegion, self.KzTransmissionRegion,
                self.layerStack, self.ExyzReflected, self.ExyzTransmitted);
        assert_almost_equal(RActual, RCalculated, self.absoluteTolerance, self.relativeTolerance);
        assert_almost_equal(TActual, TCalculated, self.absoluteTolerance, self.relativeTolerance);

    def testPMatrix(self):
        PActual = complexArray([
            [0.212504, 0.499373],
            [-0.909798, -0.212504]])
        layer = self.layerStack.internal_layers[0]
        PCalculated = layer.P_matrix()
        assert_almost_equal(PActual, PCalculated, self.absoluteTolerance, self.relativeTolerance,
                "P matrix layer 1");

    def testQMatrix(self):
        QActual = complexArray([[0.4250, 0.9987],[-1.8196, -0.4250]]);
        layer = self.layerStack.internal_layers[0]
        QCalculated = layer.Q_matrix()
        assert_almost_equal(QActual, QCalculated, self.absoluteTolerance, self.relativeTolerance,
                "Q matrix Layer 1");

        QActual = complexArray([[0.1417, 0.6662],[-0.9399, -0.1417]]);
        layer = self.layerStack.internal_layers[1]
        QCalculated = layer.Q_matrix()
        assert_almost_equal(QActual, QCalculated, self.absoluteTolerance, self.relativeTolerance,
                "Q matrix layer 2")

    def testLambdaMatrix(self):
        OActual = complexArray([[0 + 0.9046j, 0+0j],[0+0j,0+0.9046j]]);
        layer = self.layerStack.internal_layers[0]
        OCalculated = layer.lambda_matrix();
        assert_almost_equal(OActual, OCalculated, self.absoluteTolerance, self.relativeTolerance);

        OActual = complexArray([[0 + 1.3485j, 0+0j],[0+0j,0+1.3485j]]);
        layer = self.layerStack.internal_layers[1]
        OCalculated = layer.lambda_matrix();
        assert_almost_equal(OActual, OCalculated, self.absoluteTolerance, self.relativeTolerance);

    def testVMatrix(self):
        layer = self.layerStack.gapLayer
        (VCalculated, W, _, X) = layer.VWLX_matrices()
        VActual = complexArray([[0 - 0.4250j, 0 - 1.1804j], [0 + 2.0013j, 0 + 0.4250j]]);
        assert_almost_equal(VActual, VCalculated, self.absoluteTolerance, self.relativeTolerance);

        layer = self.layerStack.internal_layers[0]
        (VCalculated, W, _, X) = layer.VWLX_matrices()
        VActual = complexArray([[0-0.4698j,0-1.1040j],[0+2.0114j,0+0.4698j]]);
        assert_almost_equal(VActual, VCalculated, self.absoluteTolerance, self.relativeTolerance);

        layer = self.layerStack.internal_layers[1]
        (VCalculated, W, _, X) = layer.VWLX_matrices()
        VActual = complexArray([[0-0.1051j,0-0.4941j],[0+0.6970j,0+0.1051j]]);
        assert_almost_equal(VActual, VCalculated, self.absoluteTolerance, self.relativeTolerance);

        layer = self.layerStack.incident_layer
        (VCalculated, W_ref, _, X) = layer.VWLX_matrices()
        VActual = complexArray([
            [0 - 0.5017j, 0 - 0.8012j],
            [0 + 1.7702j, 0 + 0.5017j]]);
        assert_almost_equal(VActual, VCalculated, self.absoluteTolerance, self.relativeTolerance);

    def testXMatrix(self):
        layer = self.layerStack.internal_layers[0]
        (V, W, _, XCalculated) = layer.VWLX_matrices()
        XActual = complexArray([[0.1493+0.9888j, 0+0j],[0+0j,0.1493+0.9888j]]);
        assert_almost_equal(XActual, XCalculated, self.absoluteTolerance, self.relativeTolerance);

        layer = self.layerStack.internal_layers[1]
        (V, W, _, XCalculated) = layer.VWLX_matrices()
        XActual = complexArray([[-0.4583 - 0.8888j, 0+0j],[0+0j, -0.4583 - 0.8888j]]);
        assert_almost_equal(XActual, XCalculated, self.absoluteTolerance, self.relativeTolerance);


    def testAMatrix(self):
        W1 = np.identity(2);
        Wg = np.identity(2);
        V1 = complexArray([[0 - 0.4698j, 0 - 1.1040j],[0 + 2.0114j, 0 + 0.4698j]]);
        Vg = complexArray([[0 - 0.4250j, 0 - 1.1804j], [0 + 2.0013j, 0 + 0.4250j]]);

        ACalculated = A_matrix(W1, Wg, V1, Vg);
        AActual = self.ALayer1
        assert_almost_equal(AActual, ACalculated, self.absoluteTolerance, self.relativeTolerance);

        W2 = complexIdentity(2);
        Wg = complexIdentity(2);
        V2 = complexArray([[0 - 0.1051j, 0 - 0.4941j],[0 + 0.6970j, 0 + 0.1051j]]);
        Vg = complexArray([[0 - 0.4250j, 0 - 1.1804j],[0 + 2.0013j, 0 + 0.4250j]]);

        ACalculated = A_matrix(W2, Wg, V2, Vg);
        AActual = self.ALayer2
        assert_almost_equal(AActual, ACalculated, self.absoluteTolerance, self.relativeTolerance);

    def testBMatrix(self):
        W1 = complexIdentity(2);
        Wg = complexIdentity(2);
        V1 = complexArray([[0 - 0.4698j, 0 - 1.1040j],[0 + 2.0114j, 0 + 0.4698j]]);
        Vg = complexArray([[0 - 0.4250j, 0 - 1.1804j], [0 + 2.0013j, 0 + 0.4250j]]);
        BCalculated = B_matrix(W1, Wg, V1, Vg);
        BActual = complexArray([[-0.0049, 0.0427],[0.0427, -0.0873]]);
        assert_almost_equal(BActual, BCalculated, self.absoluteTolerance, self.relativeTolerance);

        W2 = complexIdentity(2);
        Wg = complexIdentity(2);
        V2 = complexArray([[0 - 0.1051j, 0 - 0.4941j],[0 + 0.6970j, 0 + 0.1051j]]);
        Vg = complexArray([[0 - 0.4250j, 0 - 1.1804j],[0 + 2.0013j, 0 + 0.4250j]]);

        BCalculated = B_matrix(W2, Wg, V2, Vg);
        BActual = complexArray([[-1.8324, -0.2579],[-0.2579, -1.3342]]);
        assert_almost_equal(BActual, BCalculated, self.absoluteTolerance, self.relativeTolerance);

    def testDiMatrix(self):
        absoluteTolerance = 0.003;# D matrix has some very small values after multiple matrix mult.
        relativeTolerance = 0.1; # relative error is huge on 5e-4 values. Overriding.

        A = complexArray([[2.0049, -0.0427], [-0.0427, 2.0873]]);
        B = complexArray([[-0.0049, 0.0427], [0.0427, -0.0873]]);
        X = complexArray([[0.1493 + 0.9888j, 0+0j],[0+0j, 0.4193 + 0.9888j]]);
        DCalculated = D_matrix(A, B, X);
        DActual = complexArray([[2.0057 - 0.0003j, -0.0445 + 0.0006j],[-0.0445 + 0.0006j, 2.0916 - 0.0013j]])
        assert_almost_equal(DActual, DCalculated, absoluteTolerance, relativeTolerance);

        # LAYER 2 DATA
        # Since now we have the d-matrix to higher precision we can test it more strongly.
        A = complexArray([[3.8324, 0.2579],[0.2579, 3.3342]]);
        B = complexArray([[-1.8324, -0.2579], [-0.2579, -1.3342]]);
        X = complexArray([[-0.4583 - 0.8888j, 0+0j],[0+0j, -0.4583 - 0.8888j]]);

        DCalculated = D_matrix(A, B, X);
        DActual = complexArray([[4.3436 - 0.7182j, 0.3604 - 0.1440j], [0.3604 - 0.1440j, 3.6475 - 0.4401j]]);
        assert_almost_equal(DActual, DCalculated, self.absoluteTolerance, self.relativeTolerance);

    def testScatteringMatrixFromRaw(self):
        SMatrixLayer1Calculated = calculateInternalSMatrixFromRaw(self.ALayer1, self.BLayer1,
                self.XLayer1, self.DLayer1)
        SMatrixLayer2Calculated = calculateInternalSMatrixFromRaw(self.ALayer2, self.BLayer2,
                self.XLayer2, self.DLayer2)

        S11Actual = self.S11Layer1
        S11Calculated = SMatrixLayer1Calculated[0,0];
        assert_almost_equal(S11Actual, S11Calculated, self.absoluteTolerance, self.relativeTolerance,
                "S11 for Layer 1");
        S11Calculated = SMatrixLayer2Calculated[0,0];
        S11Actual = self.S11Layer2
        assert_almost_equal(S11Actual, S11Calculated, self.absoluteTolerance, self.relativeTolerance,
                "S11 for Layer 2");

        S12Actual = self.S12Layer1
        S12Calculated = SMatrixLayer1Calculated[0,1];
        assert_almost_equal(S12Actual, S12Calculated, self.absoluteTolerance, self.relativeTolerance,
                "S12 for Layer 1");
        S12Actual = self.S12Layer2
        S12Calculated = SMatrixLayer2Calculated[0,1];
        assert_almost_equal(S12Actual, S12Calculated, self.absoluteTolerance, self.relativeTolerance,
                "S12 for Layer 2");

        S21Actual = self.S21Layer1
        S21Calculated = SMatrixLayer1Calculated[1,0];
        assert_almost_equal(S21Actual, S21Calculated, self.absoluteTolerance, self.relativeTolerance,
                "S21 for Layer 1");
        S21Actual = self.S21Layer2
        S21Calculated = SMatrixLayer2Calculated[1,0];
        assert_almost_equal(S21Actual, S21Calculated, self.absoluteTolerance, self.relativeTolerance,
                "S21 for Layer 2");

        S22Actual = self.S22Layer1
        S22Calculated = SMatrixLayer1Calculated[1,1];
        assert_almost_equal(S22Actual, S22Calculated, self.absoluteTolerance, self.relativeTolerance,
                "S22 for Layer 1");
        S22Actual = self.S22Layer2
        S22Calculated = SMatrixLayer2Calculated[1,1];
        assert_almost_equal(S22Actual, S22Calculated, self.absoluteTolerance, self.relativeTolerance,
                "S22 for Layer 2");

    def testDRedhefferMatrix(self):
        SA = self.transparentSMatrix
        SB = self.SLayer1
        DRedhefferMatrixActual = complexArray([[1,0],[0,1]])
        DRedhefferMatrixCalculated = D_matrix_redheffer(SA, SB)
        assert_almost_equal(DRedhefferMatrixActual, DRedhefferMatrixCalculated, self.absoluteTolerance,
                            self.relativeTolerance, "Layer 1 D matrix")

        SA = self.SLayer1
        SB = self.SLayer2
        DRedhefferMatrixActual = complexArray([
            [0.1506 + 0.9886j, -0.0163 - 0.0190j],
            [-0.0163 - 0.0190j, 0.1822 + 1.0253j]]);
        DRedhefferMatrixCalculated = D_matrix_redheffer(SA, SB)
        assert_almost_equal(DRedhefferMatrixActual, DRedhefferMatrixCalculated, self.absoluteTolerance,
                            self.relativeTolerance, "Layer 2 D matrix")

    def testFRedhefferMatrix(self):
        SA = self.transparentSMatrix
        SB = self.SLayer1
        FRedhefferMatrixActual = complexArray([
            [0.1490 + 0.9880j, 0.0005 + 0.0017j],
            [0.0005 + 0.0017j, 0.148 + 0.9848j]]);
        FRedhefferMatrixCalculated = F_matrix(SA, SB)
        assert_almost_equal(FRedhefferMatrixActual, FRedhefferMatrixCalculated, self.absoluteTolerance,
                            self.relativeTolerance, "Layer 1 F matrix")

        SA = self.SLayer1
        SB = self.SLayer2
        FRedhefferMatrixActual = complexArray([
            [-0.2117 - 0.6413j, 0.0471 + 0.0518j],
            [0.0471 + 0.0518j, -0.3027 - 0.7414j]]);
        FRedhefferMatrixCalculated = F_matrix(SA, SB)
        assert_almost_equal(FRedhefferMatrixActual, FRedhefferMatrixCalculated, self.absoluteTolerance,
                            self.relativeTolerance, "Layer 2 F matrix")

    def testRedhefferProduct(self):
        SA = self.transparentSMatrix
        SB = self.SLayer1
        SABActual = self.SLayer1
        SABCalculated = redheffer_product(SA, SB)
        assert_almost_equal(SABActual, SABCalculated, self.absoluteTolerance, self.relativeTolerance,
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
        SABCalculated = redheffer_product(SA, SB)
        assert_almost_equal(SABActual, SABCalculated, self.absoluteTolerance, self.relativeTolerance,
                "Redheffer product with Layer 1 Layer 2")

    def testSMatrixFromFundamentals(self):
        SiActual = self.SLayer1
        layer = self.layerStack.internal_layers[0]
        SiCalculated = layer.S_matrix()
        assert_almost_equal(SiActual, SiCalculated,
                            self.absoluteTolerance, self.relativeTolerance, "S Matrix layer 1")

        SiActual = self.SLayer2
        layer = self.layerStack.internal_layers[1]
        SiCalculated = layer.S_matrix()
        assert_almost_equal(SiActual, SiCalculated,
                            self.absoluteTolerance, self.relativeTolerance, "S Matrix layer 1")

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
        assert_almost_equal(SReflectionRegionActual, SReflectionRegionCalculated,
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
        assert_almost_equal(STransmissionRegionActual, STransmissionRegionCalculated,
                            absoluteTolerance, relativeTolerance, "S Matrix layer 1")

    def testReflectionRegionSMatrixFromFundamentals(self):
        SActual = self.SReflectionRegion
        layer = self.layerStack.incident_layer
        SCalculated = layer.S_matrix()
        assert_almost_equal(SActual, SCalculated,
                            self.absoluteTolerance, self.relativeTolerance, "S Matrix layer 1")

    def testTransmissionRegionSMatrixFromFundamentals(self):
        SActual = self.STransmissionRegion
        layer = self.layerStack.transmission_layer
        SCalculated = layer.S_matrix()
        assert_almost_equal(SActual, SCalculated,
                            self.absoluteTolerance, self.relativeTolerance, "S Matrix layer 1")

    def testCalculateTEMReflectionCoefficientsFromXYZ(self):
        #rTEMActual = complexArray([-0.418308+0.183386j, -0.222488 - 0.426831j])
        rTEMActual = complexArray([-0.591577 + 0.259348j, -0.60363 + 0.314646j])
        rTEMCalculated = calculateTEMReflectionCoefficientsFromXYZ(self.source, self.rx, self.ry, self.rz)
        assert_almost_equal(rTEMActual, rTEMCalculated,
                            self.absoluteTolerance, self.relativeTolerance, "TEM coefficients")

        # Next we test a super-simple interface against our analytic equations.

    def test_set_gap_W_matrix(self):
        self.layerStack.set_gap_layer()
        W_actual = self.layerStack.Wg
        W_desired = self.WGap
        assert_almost_equal(W_actual, W_desired,
                            self.absoluteTolerance, self.relativeTolerance)

    def test_set_gap_V_matrix(self):
        self.layerStack.set_gap_layer()
        V_actual = self.layerStack.Vg
        V_desired = self.VGap
        assert_almost_equal(V_actual, V_desired,
                            self.absoluteTolerance, self.relativeTolerance)

    @classmethod
    def setUpClass(self):
        deg = pi / 180
        self.absoluteTolerance = 1e-4
        self.relativeTolerance = 1e-3
        self.theta = 57 * deg
        self.phi = 23 * deg
        self.pTE = 1
        self.pTM = 1j
        self.wavelength = 2.7
        erReflectionRegion = 1.4
        urReflectionRegion = 1.2
        erTransmissionRegion = 1.8
        urTransmissionRegion = 1.6
        erLayer1 = 2.0
        urLayer1 = 1.0
        erLayer2 = 1.0
        urLayer2 = 3.0

        reflectionLayer = Layer(erReflectionRegion, urReflectionRegion)
        self.source = Source(wavelength=self.wavelength, theta=self.theta, phi=self.phi,
                pTEM=[self.pTE, self.pTM], layer=reflectionLayer)
        thicknessLayer1 = 0.25*self.source.wavelength
        thicknessLayer2 = 0.5*self.source.wavelength

        transmissionLayer = Layer(erTransmissionRegion, urTransmissionRegion)
        layer1 = Layer(er=erLayer1, ur=urLayer1, thickness=thicknessLayer1)
        layer2 = Layer(er=erLayer2, ur=urLayer2, thickness=thicknessLayer2)
        self.layerStack = LayerStack(layer1, layer2, incident_layer=reflectionLayer, transmission_layer=transmissionLayer)
        self.layerStack.source = self.source

        self.Kx = reflectionLayer.n* sin(self.theta) * cos(self.phi)
        self.Ky = reflectionLayer.n* sin(self.theta) * sin(self.phi)
        self.layerStack.Kx = self.Kx
        self.layerStack.Ky = self.Ky
        self.layerStack.set_gap_layer()
        self.KzReflectionRegion = 0.705995
        self.KzTransmissionRegion = 1.3032
        self.KzLayer1 = 0.9046
        self.KzLayer2 = 1.3485

        self.ExReflected = 0.0519 - 0.2856j
        self.EyReflected = -0.4324 + 0.0780j
        self.EzReflected = 0.1866 + 0.3580j
        self.ExyzReflected = complexArray([self.ExReflected, self.EyReflected, self.EzReflected])
        self.ExTransmitted = -0.0101 + 0.3577j
        self.EyTransmitted = 0.4358 - 0.0820j
        self.EzTransmitted = -0.1343 - 0.2480j
        self.ExyzTransmitted = complexArray([self.ExTransmitted, self.EyTransmitted, self.EzTransmitted])

        self.rx = 0.0519 - 0.2856j
        self.ry = -0.4324 + 0.0780j
        self.rz = 0.1866 + 0.3580j
        self.tx = -0.0101 + 0.3577j
        self.ty = 0.4358 - 0.0820j
        self.tz = -0.1343 - 0.2480j
        self.R = 0.4403
        self.T = 0.5597

        self.KzGap = 1
        self.QGap = complexArray([
            [0.42500709550286847, -0.001253971432314538],
            [-0.819595191248642, -0.42500709550286847]
        ])
        self.WGap = complexIdentity(2)
        self.VGap = complexArray([
            [0 - 0.4250j, 0 - 1.1804j],
            [0 + 2.0013j, 0 + 0.4250j]])

        self.ALayer1 = complexArray([[2.0049, -0.0427], [-0.0427, 2.0873]]);
        self.BLayer1 = complexArray([[-0.0049, 0.0427], [0.0427, -0.0873]]);
        self.XLayer1 = complexArray([[0.1493 + 0.9888j, 0+0j],[0+0j, 0.1493 + 0.9888j]]);
        self.DLayer1 = complexArray([
            [2.0057 - 0.0003j, -0.0445 + 0.0006j],
            [-0.0445 + 0.0006j, 2.0916 - 0.0013j]]);
        self.ALayer2 = complexArray([[3.8324, 0.2579],[0.2579, 3.3342]]);
        self.BLayer2 = complexArray([[-1.8324, -0.2579], [-0.2579, -1.3342]]);
        self.XLayer2 = complexArray([[-0.4583 - 0.8888j, 0+0j],[0+0j, -0.4583 - 0.8888j]]);
        self.DLayer2 = complexArray([
            [4.3436 - 0.7182j, 0.3604 - 0.1440j],
            [0.3604 - 0.1440j, 3.6475 - 0.4401j]]);
        self.S11Layer1 =  complexArray([
            [0.0039 - 0.0006j, -0.0398 + 0.0060j],
            [-0.0398 + 0.0060j, 0.0808 - 0.0121j]])
        self.S11Layer2 = complexArray([
            [0.6997 - 0.2262j, 0.0517 - 0.0014j],
            [0.0517-0.0014j, 0.5998 - 0.2235j]]);
        self.S12Layer1 = complexArray([
            [0.1490 + 0.9880j, 0.0005 + 0.0017j],
            [0.0005 + 0.0017j, 0.1480 + 0.9848j]]);
        self.S12Layer2 = complexArray([
            [-0.2093 - 0.6406j, 0.0311 + 0.0390j],
            [0.0311 + 0.0390j, -0.2693 - 0.7160j]]);
        self.S21Layer1 = complexArray([
            [0.1490 + 0.9880j, 0.0005 + 0.0017j],
            [0.0005 + 0.0017j, 0.1480 + 0.9848j]]);
        self.S21Layer2= complexArray([
            [-0.2093 - 0.6406j, 0.0311 + 0.0390j],
            [0.0311 + 0.0390j, -0.2693 - 0.7160j]]);
        self.S22Layer1 = complexArray([
            [0.0039 - 0.0006j, -0.0398 + 0.0060j],
            [-0.0398 + 0.0060j, 0.0808 - 0.0121j]]);
        self.S22Layer2 = complexArray([
            [0.6997 - 0.2262j, 0.0517 - 0.0014j],
            [0.0517-0.0014j, 0.5998 - 0.2235j]]);
        self.SLayer1 = complexArray([[self.S11Layer1, self.S12Layer1],[self.S21Layer1, self.S22Layer1]])
        self.SLayer2 = complexArray([[self.S11Layer2, self.S12Layer2],[self.S21Layer2, self.S22Layer2]])

        self.transparentSMatrix = complexZeros((2,2,2,2));
        self.transparentSMatrix[1,0] = complexIdentity(2);
        self.transparentSMatrix[0,1] = complexIdentity(2);

        self.SReflectionRegion = complexZeros((2,2,2,2));
        self.SReflectionRegion[0,0] = complexArray([
            [-0.0800, 0.0761],
            [0.0761, -0.2269]]);
        self.SReflectionRegion[0,1] = complexArray([
            [1.0800, -0.0761],
            [-0.0761, 1.2269]]);
        self.SReflectionRegion[1,0] = complexArray([
            [0.9200, 0.0761],
            [0.0761, 0.7731]]);
        self.SReflectionRegion[1,1] = complexArray([
            [0.0800, -0.0761],
            [-0.0761, 0.2269]]);

        self.STransmissionRegion = complexZeros((2,2,2,2));
        self.STransmissionRegion[0,0] = complexArray([
            [0.2060, 0.0440],
            [0.0440, 0.1209]]);
        self.STransmissionRegion[0,1] = complexArray([
            [0.7940, -0.0440],
            [-0.0440, 0.8791]]);
        self.STransmissionRegion[1,0] = complexArray([
            [1.2060, 0.0440],
            [0.0440, 1.1209]]);
        self.STransmissionRegion[1,1] = complexArray([
            [-0.2060, -0.0440],
            [-0.0440, -0.1209]]);


if __name__ == '__main__':
    unittest.main()
