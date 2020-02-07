from shorthand import *
from matrices import *
from results import *
from harmonics import *
from shorthandTest import *
from matrixParser import *
from copy import deepcopy

class RCWASolver:
    def __init__(self, layerStack, source, numberHarmonics):
        self.numberHarmonics = numberHarmonics
        self.layerStack = layerStack
        self.layerStack.setConvolutionMatrix(numberHarmonics)
        self.source = source

        self.baseCrystalLayer = self.layerStack.extractCrystalLayer()
        self.baseCrystal = self.layerStack.internalLayer[self.baseCrystalLayer].crystal
        self.setupKMatrices()
        self.KDimension = self.Kx.shape[0]
        self.scatteringElementDimension = self.KDimension * 2
        self.scatteringElementShape = (self.scatteringElementDimension, self.scatteringElementDimension)
        self.layerStack.setGapLayer(self.source.kIncident[1], self.source.kIncident[1])
        self.setupGapMatrices()
        self.setupReflectionTransmissionMatrices()

        self.ClearSolution()
        self.Si = [None for _ in range(len(self.layerStack.internalLayer))]
        self.results = []

    def setupKMatrices(self):
        self.Kx = generateKxMatrix(self.source.kIncident, self.baseCrystal, self.numberHarmonics)
        self.Ky = generateKyMatrix(self.source.kIncident, self.baseCrystal, self.numberHarmonics)
        self.KzReflectionRegion = calculateKzBackward(self.Kx, self.Ky, self.layerStack.reflectionLayer)
        self.KzTransmissionRegion = calculateKzForward(self.Kx, self.Ky, self.layerStack.transmissionLayer)
        self.KzGapRegion = calculateKzForward(self.Kx, self.Ky, self.layerStack.gapLayer)

    def Solve(self):
        self.ClearSolution()
        self.calculateDeviceSMatrix()
        self.calculateGlobalSMatrix()
        self.calculateRTQuantities()
        self.packageResults()

    def calculateRTQuantities(self):
        self.rx, self.ry, self.rz = calculateReflectionCoefficient(self.SGlobal, self.Kx, self.Ky,
                self.KzReflectionRegion, self.WReflectionRegion, self.source, self.numberHarmonics)
        self.tx, self.ty, self.tz = calculateTransmissionCoefficient(self.SGlobal, self.Kx, self.Ky,
                self.KzTransmissionRegion, self.WTransmissionRegion, self.source, self.numberHarmonics)
        self.R = calculateDiffractionReflectionEfficiency(self.rx, self.ry, self.rz, self.source,
                self.KzReflectionRegion, self.layerStack)
        self.T = calculateDiffractionTransmissionEfficiency(self.tx, self.ty, self.tz, self.source,
                self.KzTransmissionRegion, self.layerStack)
        self.RTot = np.sum(self.R)
        self.TTot = np.sum(self.T)
        self.conservation = self.RTot + self.TTot

    def packageResults(self):
        tempResults = Results()
        tempResults.rx, tempResults.ry, tempResults.rz = deepcopy((self.rx, self.ry, self.rz))
        tempResults.tx, tempResults.ty, tempResults.tz = deepcopy((self.tx, self.ty, self.tz))
        tempResults.R, tempResults.T = deepcopy((self.R, self.T))
        tempResults.RTot, tempResults.TTot, tempResults.conservation = \
                deepcopy((self.RTot, self.TTot, self.conservation))
        tempResults.crystal = deepcopy(self.baseCrystal)
        tempResults.source = deepcopy(self.source)
        tempResults.SGlobal = deepcopy(self.SGlobal)
        self.results.append(tempResults)

    def setupReflectionTransmissionMatrices(self):
        self.WReflectionRegion = complexIdentity(self.scatteringElementDimension)
        self.WTransmissionRegion = complexIdentity(self.scatteringElementDimension)

    def setupGapMatrices(self):
        self.WGap = complexIdentity(self.Kx.shape[0]*2)
        QGap = calculateQMatrix(self.Kx, self.Ky, self.layerStack.gapLayer)
        LambdaGap = calculateLambdaMatrix(self.KzGapRegion)
        self.VGap = QGap @ inv(LambdaGap)

        KzFreeSpace = numpyArrayFromFile(
                "test/matrixDataOblique/freeSpace/KzFreeSpace.txt")
        QFreeSpace = numpyArrayFromFile("test/matrixDataOblique/freeSpace/QFreeSpace.txt")
        WFreeSpace = complexIdentity(18)
        LambdaFreeSpace = numpyArrayFromSeparatedColumnsFile("test/matrixDataOblique/freeSpace/LambdaFreeSpace.txt")
        VFreeSpace = numpyArrayFromSeparatedColumnsFile("test/matrixDataOblique/freeSpace/VFreeSpace.txt")
        self.WGap = WFreeSpace # HACK HACK SHOULD BE REMOVED.
        self.VGap = VFreeSpace

    def calculateDeviceSMatrix(self):
        for i in range(len(self.layerStack.internalLayer)):
            self.Si[i] = calculateInternalSMatrix(self.Kx, self.Ky, self.layerStack.internalLayer[i],
                    self.source, self.WGap, self.VGap)
            self.SGlobal = calculateRedhefferProduct(self.SGlobal, self.Si[i])

    def calculateGlobalSMatrix(self):
        self.STransmission = calculateTransmissionRegionSMatrix(self.Kx, self.Ky, self.layerStack,
                self.WGap, self.VGap)
        self.SReflection = calculateReflectionRegionSMatrix(self.Kx, self.Ky, self.layerStack,
                self.WGap, self.VGap)
        self.SGlobal = calculateRedhefferProduct(self.SGlobal, self.STransmission)
        self.SGlobal = calculateRedhefferProduct(self.SReflection, self.SGlobal)

    def ClearSolution(self):
        self.SGlobal = generateTransparentSMatrix(self.scatteringElementShape)
        self.rx, self.ry, self.rz = None, None, None
        self.tx, self.ty, self.tz = None, None, None
        self.R, self.T, self.RTot, self.TTot, self.CTot = None, None, None, None, None
