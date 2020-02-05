from shorthand import *
from convolution import *

# WHEN WE GET BACK - THIS IS WHERE WE SHOULD PUT OUR CONVOLUTION MATRIX.
class Layer:
    def __init__(self, er=1, ur=1, L=0, crystal=None, numberHarmonics=None):
        self.er = er
        self.ur = ur
        self.L = L
        self.n = sqrt(er*ur)
        self.crystal = crystal
        if (crystal is not None) and (numberHarmonics is not None):
            self.setConvolutionMatrices(numberHarmonics)

    def setConvolutionMatrices(self, numberHarmonics):
        self.er = self.generateConvolutionMatrix(self.crystal.permittivityCellData, numberHarmonics)
        self.ur = self.generateConvolutionMatrix(self.crystal.permeabilityCellData, numberHarmonics)

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
        zeroHarmonicsLocation = [math.floor(Nx/2), math.floor(Ny/2), math.floor(Nz/2)]

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
                                desiredHarmonics = [prow - pcol, qrow - qcol, rrow - rcol]

                                # Get those harmonic locations from the zero harmonic location.
                                desiredHarmonicsLocation = zeroHarmonicsLocation + desiredHarmonics

                                convolutionMatrix[row][col] = \
                                    cellFourierRepresentation[desiredHarmonicsLocation[0]][desiredHarmonicsLocation[1]][desiredHarmonicsLocation[2]];
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

class Source:
    def __init__(self, wavelength=1,theta=0, phi=0, pTEM=[1,0], layer=freeSpaceLayer):
        self.wavelength = wavelength
        self.k0 = 2*pi / wavelength
        self.theta = theta
        self.phi = phi
        self.pTE = pTEM[0]
        self.pTM = pTEM[1]
        self.setKIncident(layer)
        self.setTEMVectors()
        self.pX = (self.pTE*self.aTE + self.pTM*self.aTM)[0]
        self.pY = (self.pTE*self.aTE + self.pTM*self.aTM)[1]

    def setTEMVectors(self):
        deviceNormalUnitVector = complexArray([0,0,-1]);
        epsilon = 1e-3;

        if(abs(self.kIncident[0]) < epsilon and abs(self.kIncident[1]) < epsilon):
            self.aTE = np.array([0,1,0]);
        else:
            self.aTE = - np.cross(deviceNormalUnitVector, self.kIncident);
            self.aTE = self.aTE / norm(self.aTE);

        self.aTM = (np.cross(self.aTE, self.kIncident));
        self.aTM = (self.aTM / norm(self.aTM));

    def setKIncident(self, layer):
        n = sqrt(layer.er*layer.ur);
        kx = n * sin(self.theta) * cos(self.phi);
        ky = n * sin(self.theta) * sin(self.phi);
        kz = n * cos(self.theta);
        self.kIncident = complexArray([kx, ky, kz])

zeroSource = Source(float("inf"))
