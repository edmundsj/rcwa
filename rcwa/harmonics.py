from rcwa.shorthand import *

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

def generateKxMatrix(source, crystal, numberHarmonics):
    if crystal is not None:
        if crystal.dimensions == 2:
            return generateKxMatrix2D(source, crystal, numberHarmonics[0:2])
        else:
            raise NotImplementedError
    else:
        return source.kIncident[0]

def generateKxMatrix2D(source, crystal, numberHarmonics):
    matrixSize = np.prod(numberHarmonics)
    matrixShape = (matrixSize, matrixSize);
    KxMatrix = complexZeros(matrixShape)

    (T1, T2) = np.array(crystal.reciprocalLatticeVectors) / source.k0
    (incidentWaveVectorx, T1x, T2x) = getXComponents(source.kIncident, T1, T2);
    (minHarmonicT1, minHarmonicT2) = calculateMinHarmonic(numberHarmonics)
    (maxHarmonicT1, maxHarmonicT2) = calculateMaxHarmonic(numberHarmonics)

    diagonalIndex = 0;
    for desiredHarmonicT2 in range(minHarmonicT2, maxHarmonicT2 + 1):
        for desiredHarmonicT1 in range(minHarmonicT1, maxHarmonicT1 + 1):

            KxMatrix[diagonalIndex][diagonalIndex] = incidentWaveVectorx - \
                    desiredHarmonicT1*T1x - desiredHarmonicT2*T2x
            diagonalIndex += 1;

    return KxMatrix

def generateKyMatrix(source, crystal, numberHarmonics):
    if crystal is not None:
        if crystal.dimensions == 2:
            return generateKyMatrix2D(source, crystal, numberHarmonics[0:2])
        else:
            raise NotImplementedError
    else:
        return source.kIncident[1]

def generateKyMatrix2D(source, crystal, numberHarmonics):
    matrixSize = np.prod(numberHarmonics)
    matrixShape = (matrixSize, matrixSize);
    KyMatrix = complexZeros(matrixShape)

    (T1, T2) = np.array(crystal.reciprocalLatticeVectors) / source.k0
    (incidentWaveVectory, T1y, T2y) = getYComponents(source.kIncident, T1, T2);
    (minHarmonicT1, minHarmonicT2) = calculateMinHarmonic(numberHarmonics)
    (maxHarmonicT1, maxHarmonicT2) = calculateMaxHarmonic(numberHarmonics)

    diagonalIndex = 0;
    for desiredHarmonicT2 in range(minHarmonicT2, maxHarmonicT2 + 1):
        for desiredHarmonicT1 in range(minHarmonicT1, maxHarmonicT1 + 1):

            KyMatrix[diagonalIndex][diagonalIndex] = incidentWaveVectory - \
                    desiredHarmonicT1*T1y - desiredHarmonicT2*T2y
            diagonalIndex += 1;

    return KyMatrix;
