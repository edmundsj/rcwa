from shorthand import *

def calculateZeroHarmonicLocation(numberHarmonics):
    zeroHarmonicLocations = [];
    for num in numberHarmonics:
        zeroHarmonicLocations.append(math.floor(num / 2));

    return zeroHarmonicLocations;

def calculateMinHarmonic(numberHarmonics):
    minHarmonics = [];
    for num in numberHarmonics:
        minHarmonics.append(- math.floor(num / 2));

    return minHarmonics;

def calculateMaxHarmonic(numberHarmonics):
    maxHarmonics = [];
    for num in numberHarmonics:
        if(num % 2 == 0):
            maxHarmonics.append(math.floor(num / 2) - 1);
        else:
            maxHarmonics.append(math.floor(num / 2));

    return maxHarmonics;


def getXComponents(*args):
    xComponents = [];
    for a in args:
        if(a.shape == (3,) or a.shape == (2,)): # element is a row vector
            xComponents.append(a[0]);
        elif(a.shape == (3,1) or a.shape == (2,1)): # element is a column vector
            xComponents.append(a[0,0]);

    return xComponents;

def getYComponents(*args):
    yComponents = [];

    for a in args:
        if(a.shape == (3,) or a.shape == (2,)):
            yComponents.append(a[1]);
        elif(a.shape == (3,1) or a.shape == (2,1)):
            yComponents.append(a[1,0]);

    return yComponents;

# Unfortunately, the source we set up (which is defined by our Kx, Ky matrices or wavevector components)
# depends on the crystal, because we need to know which modes to use. It also depends on the material
# properties of the incident layer.
class Source:
    def __init__(self, wavelength, theta, phi, pTEM, numberHarmonics, crystal, incidentLayer):
        self.wavelength = wavelength
        self.k0 = 2*pi / wavelength
        self.theta = theta
        self.phi = phi
        self.pTEM = pTEM
        self.pXY = self.pTEM[0] * aTEM[0] + pTEM[1] * aTEM[1]
        self.numberHarmonics = numberHarmonics

        setKxMatrix()

    def setKxMatrix(incidentWaveVector, crystal, numberHarmonics):
        if crystal.dimensions is 2:
            self.Kx = generateKxMatrix2D(self.incidentWaveVector, crystal, numberHarmonics[0:2])
        else:
            raise NotImplementedError

def generateKxMatrix(incidentWaveVector, crystal, numberHarmonics):
    if crystal.dimensions is 2:
        KxMatrix = generateKxMatrix2D(incidentWaveVector, crystal, numberHarmonics[0:2])
        return KxMatrix
    else:
        raise NotImplementedError

def generateKxMatrix2D(incidentWaveVector, crystal, numberHarmonics):
    matrixSize = np.prod(numberHarmonics)
    matrixShape = (matrixSize, matrixSize);
    KxMatrix = complexZeros(matrixShape)

    (T1, T2) = crystal.reciprocalLatticeVectors
    (incidentWaveVectorx, T1x, T2x) = getXComponents(incidentWaveVector, T1, T2);
    (minHarmonicT1, minHarmonicT2) = calculateMinHarmonic(numberHarmonics)
    (maxHarmonicT1, maxHarmonicT2) = calculateMaxHarmonic(numberHarmonics)

    diagonalIndex = 0;
    for desiredHarmonicT2 in range(minHarmonicT2, maxHarmonicT2 + 1):
        for desiredHarmonicT1 in range(minHarmonicT1, maxHarmonicT1 + 1):

            KxMatrix[diagonalIndex][diagonalIndex] = incidentWaveVectorx - \
                    desiredHarmonicT1*T1x - desiredHarmonicT2*T2x
            diagonalIndex += 1;

    return KxMatrix

def generateKyMatrix(incidentWaveVector, crystal, numberHarmonics):
    if crystal.dimensions is 2:
        KyMatrix = generateKyMatrix2D(incidentWaveVector, crystal, numberHarmonics[0:2])
        return KyMatrix
    else:
        raise NotImplementedError

def generateKyMatrix2D(incidentWaveVector, crystal, numberHarmonics):
    matrixSize = np.prod(numberHarmonics)
    matrixShape = (matrixSize, matrixSize);
    KyMatrix = complexZeros(matrixShape)

    (T1, T2) = crystal.reciprocalLatticeVectors
    (incidentWaveVectory, T1y, T2y) = getYComponents(incidentWaveVector, T1, T2);
    (minHarmonicT1, minHarmonicT2) = calculateMinHarmonic(numberHarmonics)
    (maxHarmonicT1, maxHarmonicT2) = calculateMaxHarmonic(numberHarmonics)

    diagonalIndex = 0;
    for desiredHarmonicT2 in range(minHarmonicT2, maxHarmonicT2 + 1):
        for desiredHarmonicT1 in range(minHarmonicT1, maxHarmonicT1 + 1):

            KyMatrix[diagonalIndex][diagonalIndex] = incidentWaveVectory - \
                    desiredHarmonicT1*T1y - desiredHarmonicT2*T2y
            diagonalIndex += 1;

    return KyMatrix;
