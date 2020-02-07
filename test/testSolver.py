import sys
sys.path.append('core');
sys.path.append('test')

import unittest
from shorthandTest import *
from matrices import *
from fresnel import *
from matrixParser import *
from source import Source
from source import Layer
from solver import *
from results import *
from crystal import Crystal


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

    def testSolve(self):
        self.solver.Solve()
        SGlobalCalculated = self.solver.SGlobal
        SGlobalActual = self.SGlobal
        assertAlmostEqual(SGlobalActual, SGlobalCalculated,
                self.absoluteTolerance, self.relativeTolerance, "testSolver: SGlobal")

    def setUp(self):
        self.absoluteTolerance = 1e-4
        self.relativeTolerance = 1e-3

        devicePermittivityCellData = np.transpose(np.loadtxt('test/triangleData.csv', delimiter=','))
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

        thicknessLayer1 = 0.5
        thicknessLayer2 = 0.3

        # THIS IS A HACK. CRYSTAL LATTICE VECTORS NEED TO BE MULTIPLIED BY K0 IN THEIR DEFINITION.
        t1 = t1 * k0
        t2 = t2 * k0
        numberHarmonics = (3, 3)

        deviceCrystal = Crystal(devicePermittivityCellData, devicePermeabilityCellData, t1, t2)
        layer1 = Layer(crystal=deviceCrystal, L=thicknessLayer1, numberHarmonics=numberHarmonics)
        layer2 = Layer(er=6.0, ur=1.0, L=thicknessLayer2)
        layerStack = LayerStack(reflectionLayer, layer1, layer2, transmissionLayer)


        self.solver = RCWASolver(layerStack, source, numberHarmonics)

        self.Kx = np.diag(complexArray(
            [2.2035, 1.0607, -0.0822, 2.2035, 1.0607, -0.0822, 2.2035, 1.0607, -0.0822]))
        self.Ky = np.diag(complexArray(
            [1.9457, 1.9457, 1.9457, 0.6124, 0.6124, 0.6124, -0.7210, -0.7210, -0.7210]))
        self.KzReflectionRegion = numpyArrayFromFile(
                "test/matrixDataOblique/reflectionRegion/KzReflectionRegion.txt")
        self.KzTransmissionRegion = np.diag(complexArray(
            [0.5989, 2.0222, 2.2820, 1.9415, 2.7386, 2.9357, 1.9039, 2.7121, 2.9109]))
        self.KzGapRegion = numpyArrayFromFile(
                "test/matrixDataOblique/freeSpace/KzFreeSpace.txt")

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

        self.S11Layer1 = numpyArrayFromSeparatedColumnsFile("test/matrixDataOblique/layer1/S11Layer1.txt")
        self.S12Layer1 = numpyArrayFromSeparatedColumnsFile("test/matrixDataOblique/layer1/S12Layer1.txt")
        self.S21Layer1 = numpyArrayFromSeparatedColumnsFile("test/matrixDataOblique/layer1/S21Layer1.txt")
        self.S22Layer1 = numpyArrayFromSeparatedColumnsFile("test/matrixDataOblique/layer1/S22Layer1.txt")
        self.SLayer1 = complexArray([
            [self.S11Layer1, self.S12Layer1],
            [self.S21Layer1, self.S22Layer1]])

        self.S11Layer2 = numpyArrayFromSeparatedColumnsFile("test/matrixDataOblique/layer2/S11Layer2.txt")
        self.S12Layer2 = numpyArrayFromSeparatedColumnsFile("test/matrixDataOblique/layer2/S12Layer2.txt")
        self.S21Layer2 = numpyArrayFromSeparatedColumnsFile("test/matrixDataOblique/layer2/S21Layer2.txt")
        self.S22Layer2 = numpyArrayFromSeparatedColumnsFile("test/matrixDataOblique/layer2/S22Layer2.txt")
        self.SLayer2 = complexArray([
            [self.S11Layer2, self.S12Layer2],
            [self.S21Layer2, self.S22Layer2]])
        self.DLayer12= np.loadtxt('test/matrixDataOblique/layer2/D12.csv', dtype=np.cdouble)
        self.FLayer12= np.loadtxt('test/matrixDataOblique/layer2/F12.csv', dtype=np.cdouble)
