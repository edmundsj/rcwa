from rcwa.shorthand import *
from typing import Union
from numpy.typing import ArrayLike
from rcwa import Crystal, Source


def zero_harmonic(n_harmonics: Union[ArrayLike, int]):
    zeroHarmonicLocations = []
    for num in n_harmonics:
        zeroHarmonicLocations.append(math.floor(num / 2))

    return zeroHarmonicLocations


def min_harmonic(n_harmonics: Union[ArrayLike, int]):
    minHarmonics = []
    if np.isscalar(n_harmonics):
        minHarmonics = - math.floor(n_harmonics / 2)
    else:
        for num in n_harmonics:
            minHarmonics.append(- math.floor(num / 2))

    return minHarmonics


def max_harmonic(n_harmonics: Union[ArrayLike, int]):
    max_harmonics = []
    if np.isscalar(n_harmonics):
        if(n_harmonics % 2 == 0):
            max_harmonics = math.floor(n_harmonics / 2) - 1
        else:
            max_harmonics = math.floor(n_harmonics / 2)
    else:
        for num in n_harmonics:
            if(num % 2 == 0):
                max_harmonics.append(math.floor(num / 2) - 1)
            else:
                max_harmonics.append(math.floor(num / 2))

    return max_harmonics


def x_components(*args: ArrayLike) -> Union[ArrayLike, complex]:
    x = []
    for a in args:
        if a.shape == (3,) or a.shape == (2,): # element is a row vector
            x.append(a[0])
        elif a.shape == (3, 1) or a.shape == (2, 1): # element is a column vector
            x.append(a[0, 0])

    if len(x) == 1:
        x = x[0]

    return x


def y_components(*args: ArrayLike) -> Union[ArrayLike, complex]:
    y = []

    for a in args:
        if a.shape == (3,) or a.shape == (2,):
            y.append(a[1])
        elif a.shape == (3, 1) or a.shape == (2, 1):
            y.append(a[1, 0])

    if len(y) == 1:
        y = y[0]
    return y


def kx_matrix(source: Source, crystal: Crystal, n_harmonics: Union[ArrayLike, int]) -> ArrayLike:
    return _k_matrix(source, crystal, n_harmonics, component='x')


def ky_matrix(source: Source, crystal: Crystal, n_harmonics: Union[ArrayLike, int]) -> ArrayLike:
    return _k_matrix(source, crystal, n_harmonics, component='y')


def _k_matrix(source: Source, crystal: Crystal,
              n_harmonics: Union[ArrayLike, int], component: str) -> ArrayLike:
    if crystal is not None:
        if crystal.dimensions == 1:
            K_matrix = _k_matrix_1D(source, crystal, n_harmonics, component=component)
        elif crystal.dimensions == 2:
            K_matrix = _k_matrix_2D(source, crystal, n_harmonics[0:2], component=component)
        else:
            raise NotImplementedError

        return K_matrix
    else:
        if component == 'x':
            kIncident_component = x_components(source.k_incident)
        elif component == 'y':
            kIncident_component = y_components(source.k_incident)
        else:
            raise ValueError(f'Component can only be x or y, not {component}')

        return kIncident_component


def _k_matrix_1D(source: Source, crystal: Crystal,
                 n_harmonics: Union[ArrayLike, int], component: str) -> ArrayLike:
    matrixSize = np.prod(n_harmonics)
    matrixShape = (matrixSize, matrixSize)
    KMatrix = complexZeros(matrixShape)
    T1 = crystal.reciprocal_lattice_vectors[0]

    if component == 'x':
        (incidentWaveVectorxy, T1xy) = x_components(source.k_incident, T1)
    elif component == 'y':
        (incidentWaveVectorxy, T1xy) = y_components(source.k_incident, T1)
    else:
        raise ValueError

    minHarmonicT1 = min_harmonic(n_harmonics)
    maxHarmonicT1 = max_harmonic(n_harmonics)

    diagonalIndex = 0
    for desiredHarmonicT1 in range(minHarmonicT1, maxHarmonicT1 + 1):

        KMatrix[diagonalIndex][diagonalIndex] = incidentWaveVectorxy - \
                desiredHarmonicT1*T1xy
        diagonalIndex += 1

    return KMatrix


def _k_matrix_2D(source, crystal, numberHarmonics, component) -> ArrayLike:
    matrixSize = np.prod(numberHarmonics)
    matrixShape = (matrixSize, matrixSize)
    KMatrix = complexZeros(matrixShape)

    (T1, T2) = np.array(crystal.reciprocal_lattice_vectors) / source.k0
    if component == 'x':
        (incidentWaveVectorxy, T1xy, T2xy) = x_components(source.k_incident, T1, T2)
    elif component == 'y':
        (incidentWaveVectorxy, T1xy, T2xy) = y_components(source.k_incident, T1, T2)
    else:
        raise ValueError

    (minHarmonicT1, minHarmonicT2) = min_harmonic(numberHarmonics)
    (maxHarmonicT1, maxHarmonicT2) = max_harmonic(numberHarmonics)

    diagonalIndex = 0
    for desiredHarmonicT2 in range(minHarmonicT2, maxHarmonicT2 + 1):
        for desiredHarmonicT1 in range(minHarmonicT1, maxHarmonicT1 + 1):

            KMatrix[diagonalIndex][diagonalIndex] = incidentWaveVectorxy - \
                    desiredHarmonicT1*T1xy - desiredHarmonicT2*T2xy
            diagonalIndex += 1

    return KMatrix
