from shorthand import *

class Layer:
    def __init__(self, er=1, ur=1, L=0, crystal=None, numberHarmonics=None):
        self.er = er
        self.ur = ur
        self.L = L
        self.n = sqrt(er*ur)
        self.crystal = crystal
        if (crystal is not None) and (numberHarmonics is not None):
            self.setConvolutionMatrices(numberHarmonics)

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
        print(P, Q, R)
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
        self.gapLayer = Layer(1,1)
        if len(layers) is 0:
            self.reflectionLayer = freeSpaceLayer
            self.transmissionLayer = freeSpaceLayer
            self.internalLayer = []
        elif len(layers) is 1:
            self.reflectionLayer = layers[0]
            self.transmissionLayer = layers[0]
            self.internalLayer = []
        else:
            self.reflectionLayer = layers[0]
            self.transmissionLayer = layers[-1]
            self.internalLayer = list(layers[1:-1])

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
