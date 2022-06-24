from rcwa.shorthand import *

def calculateZeroHarmonicLocation(numberHarmonics):
    zeroHarmonicLocations = [];
    for num in numberHarmonics:
        zeroHarmonicLocations.append(math.floor(num / 2));

    return zeroHarmonicLocations;

def calculateMinHarmonic(numberHarmonics):
    minHarmonics = [];
    if np.isscalar(numberHarmonics):
        minHarmonics = - math.floor(numberHarmonics/2)
    else:
        for num in numberHarmonics:
            minHarmonics.append(- math.floor(num / 2));

    return minHarmonics;

def calculateMaxHarmonic(numberHarmonics):
    maxHarmonics = [];
    if np.isscalar(numberHarmonics):
        if(numberHarmonics % 2 == 0):
            maxHarmonics = math.floor(numberHarmonics / 2) - 1
        else:
            maxHarmonics = math.floor(numberHarmonics / 2)
    else:
        for num in numberHarmonics:
            if(num % 2 == 0):
                maxHarmonics.append(math.floor(num / 2) - 1)
            else:
                maxHarmonics.append(math.floor(num / 2))

    return maxHarmonics;

def getXComponents(*args):
    xComponents = [];
    for a in args:
        if(a.shape == (3,) or a.shape == (2,)): # element is a row vector
            xComponents.append(a[0]);
        elif(a.shape == (3,1) or a.shape == (2,1)): # element is a column vector
            xComponents.append(a[0,0]);

    if len(xComponents) == 1:
        xComponents = xComponents[0]

    return xComponents;

def getYComponents(*args):
    yComponents = [];

    for a in args:
        if(a.shape == (3,) or a.shape == (2,)):
            yComponents.append(a[1]);
        elif(a.shape == (3,1) or a.shape == (2,1)):
            yComponents.append(a[1,0]);

    if len(yComponents) == 1:
        yComponents = yComponents[0]
    return yComponents;

def generateKxMatrix(source, crystal, numberHarmonics):
    return generateKMatrix(source, crystal, numberHarmonics, component='x')

def generateKyMatrix(source, crystal, numberHarmonics):
    return generateKMatrix(source, crystal, numberHarmonics, component='y')

def generateKMatrix(source, crystal, numberHarmonics, component):
    if crystal is not None:
        if crystal.dimensions == 1:
            K_matrix = generateKMatrix1D(source, crystal, numberHarmonics, component=component)
        elif crystal.dimensions == 2:
            K_matrix = generateKMatrix2D(source, crystal, numberHarmonics[0:2], component=component)
        else:
            raise NotImplementedError

        return K_matrix
    else:
        if component == 'x':
            kIncident_component = getXComponents(source.kIncident)
        elif component == 'y':
            kIncident_component = getYComponents(source.kIncident)
        else:
            raise ValueError(f'Component can only be x or y, not {component}')

        return kIncident_component

def generateKMatrix1D(source, crystal, numberHarmonics, component):
    matrixSize = np.prod(numberHarmonics)
    matrixShape = (matrixSize, matrixSize);
    KMatrix = complexZeros(matrixShape)
    T1 = crystal.reciprocalLatticeVectors[0]

    if component == 'x':
        (incidentWaveVectorxy, T1xy) = getXComponents(source.kIncident, T1);
    elif component == 'y':
        (incidentWaveVectorxy, T1xy) = getYComponents(source.kIncident, T1);
    else:
        raise ValueError

    minHarmonicT1 = calculateMinHarmonic(numberHarmonics)
    maxHarmonicT1 = calculateMinHarmonic(numberHarmonics)

    diagonalIndex = 0;
    for desiredHarmonicT1 in range(minHarmonicT1, maxHarmonicT1 + 1):

        KMatrix[diagonalIndex][diagonalIndex] = incidentWaveVectorxy - \
                desiredHarmonicT1*T1xy
        diagonalIndex += 1;

    return KMatrix

def generateKMatrix2D(source, crystal, numberHarmonics, component):
    matrixSize = np.prod(numberHarmonics)
    matrixShape = (matrixSize, matrixSize);
    KMatrix = complexZeros(matrixShape)

    (T1, T2) = np.array(crystal.reciprocalLatticeVectors) / source.k0
    if component == 'x':
        (incidentWaveVectorxy, T1xy, T2xy) = getXComponents(source.kIncident, T1, T2);
    elif component == 'y':
        (incidentWaveVectorxy, T1xy, T2xy) = getYComponents(source.kIncident, T1, T2);
    else:
        raise ValueError

    (minHarmonicT1, minHarmonicT2) = calculateMinHarmonic(numberHarmonics)
    (maxHarmonicT1, maxHarmonicT2) = calculateMaxHarmonic(numberHarmonics)

    diagonalIndex = 0;
    for desiredHarmonicT2 in range(minHarmonicT2, maxHarmonicT2 + 1):
        for desiredHarmonicT1 in range(minHarmonicT1, maxHarmonicT1 + 1):

            KMatrix[diagonalIndex][diagonalIndex] = incidentWaveVectorxy - \
                    desiredHarmonicT1*T1xy - desiredHarmonicT2*T2xy
            diagonalIndex += 1;

    return KMatrix;
