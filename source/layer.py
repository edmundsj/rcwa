import context
from shorthand import *

class Layer:
    def __init__(self, er=1, ur=1, L=0, n=None, crystal=None, numberHarmonics=None, material=None):
        if n == None:
            self.er = er
            self.ur = ur
            self.n = sqrt(er*ur)
        else:
            self.er = sq(n)
            self.ur = 1
            self.n = n

        self.L = L
        self.crystal = crystal
        if crystal is not None:
            self.homogenous = False
            if numberHarmonics is not None:
                self.setConvolutionMatrix(numberHarmonics)
        else:
            self.homogenous = True

    def __eq__(self, other):
        if not isinstance(other, Layer):
            return NotImplemented

        return self.er == other.er and self.ur == other.ur and self.L == other.L \
                    and self.n == other.n and self.crystal == other.crystal

    def __str__(self):
        return f'Layer with\n\ter: {self.er}\n\tur: {self.ur}\n\tL: {self.L}\n\tn: {self.n}\n\tcrystal: {self.crystal}'

    def __repr__(self):
        return f'Layer with\n\ter: {self.er}\n\tur: {self.ur}\n\tL: {self.L}\n\tn: {self.n}\n\tcrystal: {self.crystal}'


    def setConvolutionMatrix(self, numberHarmonics):
        if self.crystal is not None:
            self.er = self.generateConvolutionMatrix(self.crystal.permittivityCellData, numberHarmonics)
            self.ur = self.generateConvolutionMatrix(self.crystal.permeabilityCellData, numberHarmonics)
        else:
            self.er = self.er * complexIdentity(prod(numberHarmonics))
            self.ur = self.ur * complexIdentity(prod(numberHarmonics))

    def generateConvolutionMatrix(self, cellData, numberHarmonics):
        dataDimension = len(cellData.shape);
        if(dataDimension == 2):
            numberHarmonics = (numberHarmonics + (1,))
        (P, Q, R) = numberHarmonics

        convolutionMatrixSize = P*Q*R;
        convolutionMatrixShape = (convolutionMatrixSize, convolutionMatrixSize);
        convolutionMatrix = complexZeros(convolutionMatrixShape)

        cellData = reshapeLowDimensionalData(cellData);
        (Nx, Ny, Nz) = cellData.shape;
        zeroHarmonicsLocation = np.array([math.floor(Nx/2), math.floor(Ny/2), math.floor(Nz/2)])

        cellFourierRepresentation = fftn(cellData);
        for rrow in range(R):
            for qrow in range(Q):
                for prow in range(P):
                    row = rrow*Q*P + qrow*P + prow;
                    for rcol in range(R):
                        for qcol in range(Q):
                            for pcol in range(P):
                                col = rcol*Q*P + qcol*P + pcol;
                                # Get the desired harmonics relative to the 0th-order harmonic.
                                desiredHarmonics = np.array([prow - pcol, qrow - qcol, rrow - rcol])

                                # Get those harmonic locations from the zero harmonic location.
                                desiredHarmonicsLocation = zeroHarmonicsLocation + desiredHarmonics

                                convolutionMatrix[row][col] = \
                                    cellFourierRepresentation[desiredHarmonicsLocation[0],
                                    desiredHarmonicsLocation[1], desiredHarmonicsLocation[2]];
        if convolutionMatrix.shape == (1, 1):
            convolutionMatrix = convolutionMatrix[0][0]

        return convolutionMatrix;

freeSpaceLayer = Layer(1,1)

class LayerStack:
    def __init__(self, *layers):
        if len(layers) == 1:
            if isinstance(layers[0], list):
                layers = layers[0]
        self.gapLayer = Layer(1,1)
        if len(layers) == 0:
            self.reflectionLayer = freeSpaceLayer
            self.transmissionLayer = freeSpaceLayer
            self.internalLayer = []
        elif len(layers) == 1:
            self.reflectionLayer = layers[0]
            self.transmissionLayer = layers[0]
            self.internalLayer = []
        else:
            self.reflectionLayer = layers[0]
            self.transmissionLayer = layers[-1]
            self.internalLayer = list(layers[1:-1])

    def __eq__(self, other):
        if not isinstance(other, LayerStack):
            return NotImplemented

        reflectionLayersSame = self.reflectionLayer == other.reflectionLayer
        transmissionLayersSame = self.transmissionLayer == other.transmissionLayer
        internalLayersSame = False
        if(len(self.internalLayer) == len(other.internalLayer)):
            for i in range(len(self.internalLayer)):
                if self.internalLayer[i] != other.internalLayer[i]:
                    break;
            internalLayersSame = True

        return internalLayersSame and reflectionLayersSame and transmissionLayersSame

    def __str__(self):
        topString= f'\nReflection Layer:\n\t' + str(self.reflectionLayer) + \
                f'\nTransmissionLayer:\n\t' + str(self.transmissionLayer) + \
                f'\nInternal Layer Count: {len(self.internalLayer)}\n'
        internalString = ''
        for layer in self.internalLayer:
            internalString += str(layer) + '\n'
        return topString + internalString

    def __repr__(self):
        return f'\nReflection Layer:\n\t' + str(self.reflectionLayer) + \
                f'\nTransmissionLayer:\n\t' + str(self.transmissionLayer) + \
                f'\nInternal Layer Count: {len(self.internalLayer)}\n'
        internalString = ''
        for layer in self.internalLayer:
            internalString += str(layer) + '\n'
        return topString + internalString

    def setGapLayer(self, kx, ky):
        self.gapLayer.er = 1 + sq(kx) + sq(ky)
        self.gapLayer.ur = 1

    # set all convolution matrices for all interior layers
    def setConvolutionMatrix(self, numberHarmonics):
        for layer in self.internalLayer:
            layer.setConvolutionMatrix(numberHarmonics)

    def extractCrystalLayer(self):
        for i in range(len(self.internalLayer)):
            if self.internalLayer[i].crystal is not None:
                return i
        return 0

emptyStack = LayerStack()
