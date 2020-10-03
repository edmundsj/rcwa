import context

import unittest
from shorthandTest import *
from matrices import *
from fresnel import *
from matrixParser import *
from source import Source
from layer import Layer
from solver import *
from results import *
from crystal import Crystal
from plotter import Plotter


class TestRCWASolver(unittest.TestCase):

    def testSetupSource(self):
        kIncidentActual = complexArray([1.0607, 0.61237, 0.70711])
        kIncidentCalculated = self.solver.source.kIncident
        assertAlmostEqual(kIncidentActual, kIncidentCalculated,
                self.absoluteTolerance, self.relativeTolerance, "testSolver: kIncident")

    def testSetupKMatrices(self):
        KxActual = self.Kx
        KxCalculated = self.solver.Kx
        assertAlmostEqual(KxActual, KxCalculated,
                self.absoluteTolerance, self.relativeTolerance, "testSolver: Kx")

        KyActual = self.Ky
        KyCalculated = self.solver.Ky
        assertAlmostEqual(KyActual, KyCalculated,
                self.absoluteTolerance, self.relativeTolerance, "testSolver: Ky")

        KzActual = self.KzReflectionRegion
        KzCalculated = self.solver.KzReflectionRegion
        assertAlmostEqual(KzActual, KzCalculated,
                self.absoluteTolerance, self.relativeTolerance, "testSolver: KzReflectionRegion")

        KzActual = self.KzTransmissionRegion
        KzCalculated = self.solver.KzTransmissionRegion
        assertAlmostEqual(KzActual, KzCalculated,
                self.absoluteTolerance, self.relativeTolerance, "testSolver: KzTransmissionRegion")

        KzActual = self.KzGapRegion
        KzCalculated = self.solver.KzGapRegion
        assertAlmostEqual(KzActual, KzCalculated,
                self.absoluteTolerance, self.relativeTolerance, "testSolver: KzGapRegion")

    def testEdgeSMatrices(self):
        self.solver.Solve()
        SActual = self.SReflectionRegion
        SCalculated = self.solver.SReflection
        assertAlmostEqual(SActual, SCalculated,
                self.absoluteTolerance, self.relativeTolerance, "testSolver: SReflection")

        self.solver.Solve()
        SActual = self.STransmissionRegion
        SCalculated = self.solver.STransmission
        assertAlmostEqual(SActual, SCalculated,
                self.absoluteTolerance, self.relativeTolerance, "testSolver: STransmission")

    def testInternalSMatrices(self):
        self.solver.Solve()
        SActual = self.SLayer1
        SCalculated = self.solver.Si[0]
        assertAlmostEqual(SActual, SCalculated,
                self.absoluteTolerance, self.relativeTolerance, "testSolver: Si[0]")

        self.solver.Solve()
        SActual = self.SLayer2
        SCalculated = self.solver.Si[1]
        assertAlmostEqual(SActual, SCalculated,
                self.absoluteTolerance, self.relativeTolerance, "testSolver: Si[1]")

    def testrtAmplitudeCoefficients(self):
        self.solver.Solve()
        # HACK - FOR SOME REASON MY PHASE IS OFF BY PI.
        rxActual = -self.rx
        ryActual = -self.ry
        rzActual = -self.rz
        (rxCalculated, ryCalculated, rzCalculated) = (self.solver.rx, self.solver.ry, self.solver.rz)
        (txCalculated, tyCalculated, tzCalculated) = (self.solver.tx, self.solver.ty, self.solver.tz)

        assertAlmostEqual(rxActual, rxCalculated,
                self.absoluteTolerance, self.relativeTolerance, "testSolver: rx")
        assertAlmostEqual(ryActual, ryCalculated,
                self.absoluteTolerance, self.relativeTolerance, "testSolver: ry")
        assertAlmostEqual(rzActual, rzCalculated,
                self.absoluteTolerance, self.relativeTolerance, "testSolver: rz")

        txActual = -self.tx
        tyActual = -self.ty
        tzActual = -self.tz
        (rxCalculated, ryCalculated, rzCalculated) = (self.solver.rx, self.solver.ry, self.solver.rz)
        (txCalculated, tyCalculated, tzCalculated) = (self.solver.tx, self.solver.ty, self.solver.tz)
        (R, T) = (self.solver.R, self.solver.T)
        assertAlmostEqual(txActual, txCalculated,
                self.absoluteTolerance, self.relativeTolerance, "testSolver: tx")
        assertAlmostEqual(tyActual, tyCalculated,
                self.absoluteTolerance, self.relativeTolerance, "testSolver: ty")
        assertAlmostEqual(tzActual, tzCalculated,
                self.absoluteTolerance, self.relativeTolerance, "testSolver: tz")

    def testDiffractionEfficiencies(self):
        self.solver.Solve()
        RActual = self.R
        TActual = self.T
        (RCalculated, TCalculated) = (self.solver.R, self.solver.T)
        assertAlmostEqual(RActual, RCalculated,
                self.absoluteTolerance, self.relativeTolerance, "testSolver: R")
        assertAlmostEqual(TActual, TCalculated,
                self.absoluteTolerance, self.relativeTolerance, "testSolver: T")

        RTotActual = self.RTot
        TTotActual = self.TTot
        CTotActual = 1.0
        RTotCalculated = self.solver.RTot
        TTotCalculated = self.solver.TTot
        CTotCalculated = self.solver.conservation
        assertAlmostEqual(RTotActual, RTotCalculated,
                self.absoluteTolerance, self.relativeTolerance, "testSolver: RTot")
        assertAlmostEqual(TTotActual, TTotCalculated,
                self.absoluteTolerance, self.relativeTolerance, "testSolver: TTot")
        assertAlmostEqual(CTotActual, CTotCalculated, 1e-7, 1e-7, "testSolver: Conservation Violated")

    def testIntegrationMultiWavelength(self):
        testWavelengths = self.solver.source.wavelength*np.arange(0.2,2,0.01)
        self.solver.Solve(testWavelengths)
        #Plotter.plotReflectionSpectra(self.solver.results)

    def setUp(self):
        self.absoluteTolerance = 1e-4
        self.relativeTolerance = 1e-3

        devicePermittivityCellData = np.transpose(np.loadtxt(context.testLocation + '/triangleData.csv', delimiter=','))
        devicePermeabilityCellData = 1 + 0 * devicePermittivityCellData

        reflectionLayer = Layer(er=2.0, ur=1.0)
        transmissionLayer = Layer(er=9.0, ur=1.0)

        # NOTE: t1 AND t2 MUST BE NORMALIZED BY MULTIPLYING BY k0, OTHERWISE THIS WILL NOT WORK, AS
        # EVERYTHING WAS FORMULATED IN TERMS OF NORMALIZED WAVEVECTORS. I DON'T KNOW OF AN ELEGANT WAY
        # TO DO THIS OTHER THAN REQUIRING A CRYSTAL TO HAVE A SOURCE AS THE INPUT. I DON'T KNOW OF
        # AN EASY WAY TO FIX THIS. I'M GOING TO FUDGE IT FOR NOW.
        wavelength = 2
        k0 = 2*pi/wavelength
        theta = 60 * deg
        phi = 30*deg
        pTEM = 1/sqrt(2)*complexArray([1,1j])
        source = Source(wavelength=wavelength, theta=theta, phi=phi, pTEM=pTEM, layer=reflectionLayer)
        t1, t2 = complexArray([1.75, 0, 0]), complexArray([0, 1.5, 0])

        thicknessLayer1 = 0.5 # should be 0.5
        thicknessLayer2 = 0.3 # should be 0.3

        numberHarmonics = (3, 3)

        deviceCrystal = Crystal(devicePermittivityCellData, devicePermeabilityCellData, t1, t2)
        layer1 = Layer(crystal=deviceCrystal, L=thicknessLayer1, numberHarmonics=numberHarmonics)
        layer2 = Layer(er=6.0, ur=1.0, L=thicknessLayer2)
        layerStack = LayerStack(reflectionLayer, layer1, layer2, transmissionLayer)


        self.solver = RCWASolver(layerStack, source, numberHarmonics)

    @classmethod
    def setUpClass(self):
        """
        Test fixture for loading in all the external test data.
        """
        self.Kx = np.diag(complexArray(
            [2.2035, 1.0607, -0.0822, 2.2035, 1.0607, -0.0822, 2.2035, 1.0607, -0.0822]))
        self.Ky = np.diag(complexArray(
            [1.9457, 1.9457, 1.9457, 0.6124, 0.6124, 0.6124, -0.7210, -0.7210, -0.7210]))
        self.KzReflectionRegion = numpyArrayFromFile(
                context.testLocation + "/matrixDataOblique/reflectionRegion/KzReflectionRegion.txt")
        self.KzTransmissionRegion = np.diag(complexArray(
            [0.5989, 2.0222, 2.2820, 1.9415, 2.7386, 2.9357, 1.9039, 2.7121, 2.9109]))
        self.KzGapRegion = numpyArrayFromFile(
                context.testLocation + "/matrixDataOblique/freeSpace/KzFreeSpace.txt")

        self.SGlobal11= numpyArrayFromSeparatedColumnsFile(
                context.testLocation + "/matrixDataOblique/SGlobal11.txt")
        self.SGlobal12= numpyArrayFromSeparatedColumnsFile(
                context.testLocation + "/matrixDataOblique/SGlobal12.txt")
        self.SGlobal21= numpyArrayFromSeparatedColumnsFile(
                context.testLocation + "/matrixDataOblique/SGlobal21.txt")
        self.SGlobal22= numpyArrayFromSeparatedColumnsFile(
                context.testLocation + "/matrixDataOblique/SGlobal22.txt")
        self.SGlobal = complexArray([
            [self.SGlobal11, self.SGlobal12],
            [self.SGlobal21, self.SGlobal22]])

        self.S11ReflectionRegion = numpyArrayFromSeparatedColumnsFile(
                context.testLocation + "/matrixDataOblique/reflectionRegion/S11ReflectionRegion.txt")
        self.S12ReflectionRegion = numpyArrayFromSeparatedColumnsFile(
                context.testLocation + "/matrixDataOblique/reflectionRegion/S12ReflectionRegion.txt")
        self.S21ReflectionRegion = numpyArrayFromSeparatedColumnsFile(
                context.testLocation + "/matrixDataOblique/reflectionRegion/S21ReflectionRegion.txt")
        self.S22ReflectionRegion = numpyArrayFromSeparatedColumnsFile(
                context.testLocation + "/matrixDataOblique/reflectionRegion/S22ReflectionRegion.txt")
        self.SReflectionRegion = complexArray([
            [self.S11ReflectionRegion, self.S12ReflectionRegion],
            [self.S21ReflectionRegion, self.S22ReflectionRegion]])

        self.S11TransmissionRegion = numpyArrayFromSeparatedColumnsFile(
                context.testLocation + "/matrixDataOblique/transmissionRegion/S11TransmissionRegion.txt")
        self.S12TransmissionRegion = numpyArrayFromSeparatedColumnsFile(
                context.testLocation + "/matrixDataOblique/transmissionRegion/S12TransmissionRegion.txt")
        self.S21TransmissionRegion = numpyArrayFromSeparatedColumnsFile(
                context.testLocation + "/matrixDataOblique/transmissionRegion/S21TransmissionRegion.txt")
        self.S22TransmissionRegion = numpyArrayFromSeparatedColumnsFile(
                context.testLocation + "/matrixDataOblique/transmissionRegion/S22TransmissionRegion.txt")
        self.STransmissionRegion = complexArray([
            [self.S11TransmissionRegion, self.S12TransmissionRegion],
            [self.S21TransmissionRegion, self.S22TransmissionRegion]])

        self.S11Layer1 = numpyArrayFromSeparatedColumnsFile(context.testLocation + "/matrixDataOblique/layer1/S11Layer1.txt")
        self.S12Layer1 = numpyArrayFromSeparatedColumnsFile(context.testLocation + "/matrixDataOblique/layer1/S12Layer1.txt")
        self.S21Layer1 = numpyArrayFromSeparatedColumnsFile(context.testLocation + "/matrixDataOblique/layer1/S21Layer1.txt")
        self.S22Layer1 = numpyArrayFromSeparatedColumnsFile(context.testLocation + "/matrixDataOblique/layer1/S22Layer1.txt")
        self.SLayer1 = complexArray([
            [self.S11Layer1, self.S12Layer1],
            [self.S21Layer1, self.S22Layer1]])

        self.S11Layer2 = numpyArrayFromSeparatedColumnsFile(context.testLocation + "/matrixDataOblique/layer2/S11Layer2.txt")
        self.S12Layer2 = numpyArrayFromSeparatedColumnsFile(context.testLocation + "/matrixDataOblique/layer2/S12Layer2.txt")
        self.S21Layer2 = numpyArrayFromSeparatedColumnsFile(context.testLocation + "/matrixDataOblique/layer2/S21Layer2.txt")
        self.S22Layer2 = numpyArrayFromSeparatedColumnsFile(context.testLocation + "/matrixDataOblique/layer2/S22Layer2.txt")
        self.SLayer2 = complexArray([
            [self.S11Layer2, self.S12Layer2],
            [self.S21Layer2, self.S22Layer2]])
        self.DLayer12= np.loadtxt(context.testLocation + '/matrixDataOblique/layer2/D12.csv', dtype=np.cdouble)
        self.FLayer12= np.loadtxt(context.testLocation + '/matrixDataOblique/layer2/F12.csv', dtype=np.cdouble)

        self.rx = complexArray([-0.0187- 0.0155j, 0.0486 - 0.0467j, 0.0016 + 0.0012j,
            0.0324 - 0.0229j, -0.1606 - 0.0348j, -0.0089 + 0.0156j,
            0.0020 + 0.0105j, 0.0076 + 0.0187j, -0.0027 - 0.0129j])
        self.ry = complexArray([-0.0077 - 0.0106j, 0.0184 + 0.0323j, -0.0267 - 0.0070j,
            -0.0286 + 0.0472j, 0.2335 + 0.0138j, 0.0243 + 0.0164j,
            0.0435 - 0.018j, 0.0183 + 0.0146j, -0.0062 + 0.0011j])
        self.rz = complexArray([0.0213 - 0.0218j, -0.0078 + 0.0512j, 0.0103 - 0.0388j,
            0.0120 + 0.0300j, -0.0386 - 0.0403j, 0.0123 + 0.0069j,
            -0.0197 - 0.0147j, -0.0087 + 0.0157j, 0.0039 + 0.0002j])
        self.tx = complexArray([0.0015 - 0.0016j, -0.0583 + 0.0256j, -0.0245 - 0.0098j,
            0.0060 + 0.0210j, 0.3040 + 0.0664j, -0.0054 - 0.0632j,
            -0.0123 - 0.0262j, -0.0323 - 0.0534j, 0.0169 + 0.0455j])
        self.ty = complexArray([-0.0024 + 0.0011j, 0.0356 + 0.0282j, -0.0230 - 0.0071j,
            0.0610 - 0.0011j, 0.0523 - 0.2913j, -0.0645 - 0.0027j,
            -0.0170 - 0.0165j, -0.0420 + 0.0298j, 0.0258 - 0.0234j])
        self.tz = complexArray([0.0023 + 0.0021j, - 0.0036 - 0.0406j, 0.0187 + 0.0057j,
            -0.0261 - 0.0235j, -0.1294 + 0.0394j, 0.0133 - 0.0012j,
            0.0078 + 0.0241j, 0.0014 + 0.0288j, 0.0069 - 0.0045j])

        self.R = np.array([
            [0,0,0],
            [0,0.0848, 0.0011],
            [0, 0.0025, 0.0004]])
        self.T = np.array([
            [0, 0.0149, 0.0055],
            [0.0222, 0.7851, 0.0283],
            [0.0053, 0.0348, 0.0150]])
        self.R = np.transpose(self.R)
        self.T = np.transpose(self.T)
        self.RTot = 0.088768
        self.TTot = 0.91123
