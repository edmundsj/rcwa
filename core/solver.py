from shorthand import *
from matrices import *
from results import *
from harmonics import *

class RCWASolver:
    def __init__(self, layerStack, source, numberHarmonics):
        self.numberHarmonics = numberHarmonics
        self.layerStack = layerStack
        self.layerStack.setConvolutionMatrix(numberHarmonics)
        self.source = source

        self.baseCrystalLayer = self.layerStack.extractCrystalLayer()
        self.baseCrystal = self.layerStack.internalLayer[self.baseCrystalLayer].crystal
        self.setupKMatrices()
        self.layerStack.setGapLayer(self.source.kIncident[1], self.source.kIncident[1])

        scatteringElementDimension = self.Kx.shape[0] * 2
        scatteringElementShape = (scatteringElementDimension, scatteringElementDimension)
        self.SGlobal = generateTransparentSMatrix(scatteringElementShape)

        # In order to calculate the gap matrices, we need the VWX function not to solve the 
        # eigenvalue problem when givien a homogenous media... This I will implement tomorrow.


    def setupKMatrices(self):
        self.Kx = generateKxMatrix(self.source.kIncident, self.baseCrystal, self.numberHarmonics)
        self.Ky = generateKyMatrix(self.source.kIncident, self.baseCrystal, self.numberHarmonics)
        self.KzReflectionRegion = calculateKzBackward(self.Kx, self.Ky, self.layerStack.reflectionLayer)
        self.KzTransmissionRegion = calculateKzForward(self.Kx, self.Ky, self.layerStack.transmissionLayer)

    def calculateDeviceMatrix(self):
        for layer in self.layerStack.internalLayer:
            Si = calculateInternalSMatrix(self.Kx, self.Ky, layer, self.source, self.Wg, self.Vg)
            self.SGlobal = redhefferStarProduct(self.SGlobal, Si)

    def Solve(self):
        calculateDeviceMatrix()
