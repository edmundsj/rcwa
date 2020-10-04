import context
from shorthand import *
from matrices import *
from harmonics import *
from matrixParser import *
from layer import Layer, LayerStack
from copy import deepcopy

class Solver:
    """ Main class that invokes all methods necessary to solve an RCWA/TMM simulation

    :param layerStack: layerStack: Stack of layers to simulate
    :param source: Source object which includes wavelength and direction information
    :param numberHarmonics: The number of harmonics in x, y to simulate (number of Fourier components). Defaults to (1, 1). Should be a 2-dimensional tuple.
    """
    def __init__(self, layerStack, source, numberHarmonics=(1,1)):
        self.numberHarmonics = numberHarmonics
        self.layerStack = layerStack
        self.layerStack.setConvolutionMatrix(numberHarmonics)
        self.source = source
        self.source.layer = layerStack.reflectionLayer

        self.baseCrystalLayer = self.layerStack.extractCrystalLayer()
        self.baseCrystal = self.layerStack.internalLayer[self.baseCrystalLayer].crystal
        self.setupKMatrices()
        self.setupGapMatrices()
        self.setupReflectionTransmissionMatrices()
        self.results = []
        self.ClearSolution()

    def setupKMatrices(self):
        """
        Sets up the Kx, Ky, and Kz matrices for solving the simulation once the source, crystal, and
        number harmonics are known.
        """
        self.Kx = generateKxMatrix(self.source, self.baseCrystal, self.numberHarmonics)
        self.Ky = generateKyMatrix(self.source, self.baseCrystal, self.numberHarmonics)
        if isinstance(self.Kx, np.ndarray):
            self.KDimension = self.Kx.shape[0]
            self.TMMSimulation = False
        else:
            self.KDimension = 1
            # Ensure that Kz for the gap layer is 1
            self.layerStack.gapLayer = Layer(er=1 + sq(self.Kx) + sq(self.Ky), ur=1, L=0)
            self.TMMSimulation = True
        self.KzReflectionRegion = calculateKzBackward(self.Kx, self.Ky, self.layerStack.reflectionLayer)
        self.KzTransmissionRegion = calculateKzForward(self.Kx, self.Ky, self.layerStack.transmissionLayer)
        self.KzGapRegion = calculateKzForward(self.Kx, self.Ky, self.layerStack.gapLayer)

        self.scatteringElementDimension = self.KDimension * 2
        self.scatteringElementShape = (self.scatteringElementDimension, self.scatteringElementDimension)

    def Solve(self, wavelengths=np.array([])):
        """
        Solves the simulation or performs a simulation sweep of the desired parameters

        :param wavelengths: list of wavelengths to simulate (can be a single number)
        """
        if wavelengths.size == 0:
            wavelengths = np.array([self.source.wavelength])
        for wavelength in wavelengths:
            self.source.layer = self.layerStack.reflectionLayer
            self.source.wavelength = wavelength # Update the source wavelength and all associated things.
            self.setupKMatrices()
            self.setupGapMatrices()
            self.setupReflectionTransmissionMatrices()
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

        if self.TMMSimulation is True:
            self.rTEM = calculateTEMReflectionCoefficientsFromXYZ(self.source, self.rx, self.ry, self.rz)

    def packageResults(self):
        """
        Packages the results from the simulation into a dictionary
        """
        tempResults = {}
        tempResults['rx'], tempResults['ry'], tempResults['rz'] = deepcopy((self.rx, self.ry, self.rz))
        tempResults['tx'], tempResults['ty'], tempResults['tz'] = deepcopy((self.tx, self.ty, self.tz))
        tempResults['R'], tempResults['T'] = deepcopy((self.R, self.T))
        tempResults['RTot'], tempResults['TTot'], tempResults['conservation'] = \
                deepcopy((self.RTot, self.TTot, self.conservation))
        tempResults['crystal'] = deepcopy(self.baseCrystal)
        tempResults['source'] = deepcopy(self.source)
        tempResults['S'] = deepcopy(self.SGlobal)

        if self.TMMSimulation is True:
            tempResults['rTE'] = self.rTEM[0]
            tempResults['rTM'] = self.rTEM[1]
            rho = tempResults['rTM'] / tempResults['rTE']
            tempResults['tanPsi'] = np.abs(rho)
            tempResults['cosDelta'] = cos(np.angle(rho))
            tempResults['delta'] = np.angle(rho)

        self.results.append(tempResults)

    def setupReflectionTransmissionMatrices(self):
        self.WReflectionRegion = complexIdentity(self.scatteringElementDimension)
        self.WTransmissionRegion = complexIdentity(self.scatteringElementDimension)

    def setupGapMatrices(self):
        self.WGap = complexIdentity(self.scatteringElementDimension)
        QGap = calculateQMatrix(self.Kx, self.Ky, self.layerStack.gapLayer)
        LambdaGap = calculateLambdaMatrix(self.KzGapRegion)
        self.VGap = QGap @ inv(LambdaGap)

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
        self.Si = [None for _ in range(len(self.layerStack.internalLayer))]
