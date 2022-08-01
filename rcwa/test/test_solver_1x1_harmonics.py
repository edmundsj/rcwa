# DATA FOR THIS FILE TAKEN FROM THE TMM BENCHMARKING DATA PROVIDED BY RUMPF.
import unittest
from rcwa import Source, Layer, Crystal, Solver, LayerStack
from rcwa.shorthand import *
from rcwa.testing import *
from rcwa.matrices import *
import numpy as np

class TestSolver1x1(unittest.TestCase):

    def testSetupSource(self):
        kIncidentActual = self.kIncident
        kIncidentCalculated = self.solver.source.k_incident
        assert_almost_equal(kIncidentActual, kIncidentCalculated,
                            self.absoluteTolerance, self.relativeTolerance, "testSolver1x1: kIncident")

    def testSetupKMatrices(self):
        KxActual = self.Kx
        KxCalculated = self.solver.Kx
        assert_almost_equal(KxActual, KxCalculated,
                            self.absoluteTolerance, self.relativeTolerance, "testSolver1x1: Kx")

        KyActual = self.Ky
        KyCalculated = self.solver.Ky
        assert_almost_equal(KyActual, KyCalculated,
                            self.absoluteTolerance, self.relativeTolerance, "testSolver1x1: Ky")

        KzActual = self.KzReflectionRegion
        KzCalculated = self.solver.KzReflectionRegion
        assert_almost_equal(KzActual, KzCalculated,
                            self.absoluteTolerance, self.relativeTolerance, "testSolver1x1: KzReflectionRegion")

        KzActual = self.KzTransmissionRegion
        KzCalculated = self.solver.KzTransmissionRegion
        assert_almost_equal(KzActual, KzCalculated,
                            self.absoluteTolerance, self.relativeTolerance, "testSolver1x1: KzTransmissionRegion")

        KzActual = self.KzGapRegion
        KzCalculated = self.solver.KzGapRegion
        assert_almost_equal(KzActual, KzCalculated,
                            self.absoluteTolerance, self.relativeTolerance, "testSolver1x1: KzGapRegion")

    def testEdgeSMatrices(self):
        self.solver.solve()
        SActual = self.SReflectionRegion
        SCalculated = self.solver.SReflection
        assert_almost_equal(SActual, SCalculated,
                            self.absoluteTolerance, self.relativeTolerance, "testSolver1x1: SReflection")

        self.solver.solve()
        SActual = self.STransmissionRegion
        SCalculated = self.solver.STransmission
        assert_almost_equal(SActual, SCalculated,
                            self.absoluteTolerance, self.relativeTolerance, "testSolver1x1: STransmission")

    def testInternalSMatrices(self):
        self.solver.solve()
        SActual = self.SLayer1
        SCalculated = self.solver.Si[0]
        assert_almost_equal(SActual, SCalculated,
                            self.absoluteTolerance, self.relativeTolerance, "testSolver1x1: Si[0]")

        SActual = self.SLayer2
        SCalculated = self.solver.Si[1]
        assert_almost_equal(SActual, SCalculated,
                            self.absoluteTolerance, self.relativeTolerance, "testSolver1x1: Si[1]")


    def testrtAmplitudeCoefficients(self):
        self.solver.solve()
        rxActual = self.rx
        ryActual = self.ry
        rzActual = self.rz
        (rxCalculated, ryCalculated, rzCalculated) = (self.solver.rx, self.solver.ry, self.solver.rz)

        assert_almost_equal(rxActual, rxCalculated,
                            self.absoluteTolerance, self.relativeTolerance, "testSolver1x1: rx")
        assert_almost_equal(ryActual, ryCalculated,
                            self.absoluteTolerance, self.relativeTolerance, "testSolver1x1: ry")
        assert_almost_equal(rzActual, rzCalculated,
                            self.absoluteTolerance, self.relativeTolerance, "testSolver1x1: rz")

        txActual = self.tx
        tyActual = self.ty
        tzActual = self.tz
        (txCalculated, tyCalculated, tzCalculated) = (self.solver.tx, self.solver.ty, self.solver.tz)
        assert_almost_equal(txActual, txCalculated,
                            self.absoluteTolerance, self.relativeTolerance, "testSolver1x1: tx")
        assert_almost_equal(tyActual, tyCalculated,
                            self.absoluteTolerance, self.relativeTolerance, "testSolver1x1: ty")
        assert_almost_equal(tzActual, tzCalculated,
                            self.absoluteTolerance, self.relativeTolerance, "testSolver1x1: tz")

        rTEActual = self.rTE
        rTMActual = self.rTM
        rTECalculated = self.solver.rTEM[0]
        rTMCalculated = self.solver.rTEM[1]
        assert_almost_equal(rTEActual, rTECalculated,
                            self.absoluteTolerance, self.relativeTolerance, "testSolver1x1: rTE")
        assert_almost_equal(rTMActual, rTMCalculated,
                            self.absoluteTolerance, self.relativeTolerance, "testSolver1x1: rTM")

    def testEllipsometryData(self):
        self.solver.solve()
        tanPsiActual = self.tanPsi
        cosDeltaActual = self.cosDelta
        deltaActual = self.delta

        tanPsiCalculated = self.solver.results['tanPsi']
        deltaCalculated = self.solver.results['delta']
        cosDeltaCalculated = self.solver.results['cosDelta']

        assert_almost_equal(tanPsiActual, tanPsiCalculated,
                            self.absoluteTolerance, self.relativeTolerance, "testSolver1x1: tanPsi")
        assert_almost_equal(deltaActual, deltaCalculated,
                            self.absoluteTolerance, self.relativeTolerance, "testSolver1x1: delta")
        assert_almost_equal(cosDeltaActual, cosDeltaCalculated,
                            self.absoluteTolerance, self.relativeTolerance, "testSolver1x1: cosDelta")

    def testRT(self):
        self.solver.solve()
        RActual = self.R
        TActual = self.T
        (RCalculated, TCalculated) = (self.solver.R, self.solver.T)
        assert_almost_equal(RActual, RCalculated,
                            self.absoluteTolerance, self.relativeTolerance, "testSolver1x1: R")
        assert_almost_equal(TActual, TCalculated,
                            self.absoluteTolerance, self.relativeTolerance, "testSolver1x1: T")

        RTotActual= self.R
        TTotActual = self.T
        (RTotCalculated, TTotCalculated) = (self.solver.RTot, self.solver.TTot)
        assert_almost_equal(RTotActual, RTotCalculated,
                            self.absoluteTolerance, self.relativeTolerance, "testSolver1x1: R")
        assert_almost_equal(TTotActual, TTotCalculated,
                            self.absoluteTolerance, self.relativeTolerance, "testSolver1x1: T")

    def testPackageResults(self):
        """ UNFINISHED """
        self.solver.solve()
        RActual = self.R
        TActual = self.T

    def testIntegrationMultiWavelength(self):
        testWavelengths = self.solver.source.wavelength*np.arange(0.2,2,0.01)
        self.solver.solve(wavelength=testWavelengths)

    def setUp(self):
        self.absoluteTolerance = 1e-4
        self.relativeTolerance = 1e-3

        reflectionLayer = Layer(er=1.4, ur=1.2)
        transmissionLayer = Layer(er=1.8, ur=1.6)

        # NOTE: t1 AND t2 MUST BE NORMALIZED BY MULTIPLYING BY k0, OTHERWISE THIS WILL NOT WORK, AS
        # EVERYTHING WAS FORMULATED IN TERMS OF NORMALIZED WAVEVECTORS. I DON'T KNOW OF AN ELEGANT WAY
        # TO DO THIS OTHER THAN REQUIRING A CRYSTAL TO HAVE A SOURCE AS THE INPUT. I DON'T KNOW OF
        # AN EASY WAY TO FIX THIS. I'M GOING TO FUDGE IT FOR NOW.
        wavelength = 2.7
        k0 = 2*pi/wavelength
        theta = 57 * deg
        phi = 23*deg
        pTEM = 1/sqrt(2)*complexArray([1,1j])
        source = Source(wavelength=wavelength, theta=theta, phi=phi, pTEM=pTEM, layer=reflectionLayer)

        thicknessLayer1 = wavelength * 0.25 # should be 0.5
        thicknessLayer2 = wavelength * 0.5 # should be 0.3

        numberHarmonics = (1,1)

        layer1 = Layer(er=2.0, ur=1.0, thickness=thicknessLayer1)
        layer2 = Layer(er=1.0, ur=3.0, thickness=thicknessLayer2)
        layerStack = LayerStack(layer1, layer2, incident_layer=reflectionLayer, transmission_layer=transmissionLayer)

        self.solver = Solver(layerStack, source, numberHarmonics)

        self.Kx = 1.0006
        self.Ky = 0.4247
        self.KzReflectionRegion = 0.705995
        self.kIncident = complexArray([self.Kx, self.Ky, self.KzReflectionRegion])
        self.KzTransmissionRegion = 1.30324
        self.KzGapRegion = 1

        self.rx = 0.0519 - 0.2856j
        self.ry= -0.4324 + 0.0780j
        self.rz = 0.1866 + 0.3580j
        self.tx = -0.0101 + 0.3577j
        self.ty = 0.4358 - 0.0820j
        self.tz = -0.1343 - 0.2480j
        self.R = 0.4403
        self.T = 0.5597

        self.tanPsi = 1.0538
        self.cosDelta = 0.997733
        self.delta = -0.0673421

        #self.rTE = -0.418308 + 0.183386j
        #self.rTM = -0.222488 - 0.426831j
        self.rTE = -0.591577 + 0.259348j
        self.rTM = -0.60363 + 0.314646j

        self.KzGap = 1
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
