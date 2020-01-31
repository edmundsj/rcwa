from shorthand import *

OUTER_BLOCK_SHAPE = (2,2);
scatteringElementShape = (2,2); # The shape of our core PQ matrices.
scatteringMatrixShape = OUTER_BLOCK_SHAPE + scatteringElementShape;
scatteringElementShape = (2,2);
scatteringElementSize = scatteringElementShape[0];

def aTEMGen(kx, ky, kz):
    deviceNormalUnitVector = complexArray([0,0,-1]); # The z-normal vector (plane wave opposite direction as surface)
    kn_vec = complexArray([kx, ky, kz]);
    epsilon = 1e-3; # a small number. If both our kx and ky vectors are tiny, our cross product will not

    if(abs(kx) < epsilon and abs(ky) < epsilon):
        aTE = np.array([0,1,0]); # This is assuming W is the identity matrix.
    else:
        aTE = - np.cross(deviceNormalUnitVector, kn_vec);
        aTE = aTE / norm(aTE);

    aTM = np.cross(aTE, kn_vec);
    aTM = aTM / norm(aTM);
    return (aTE[0:2], aTM[0:2]);

def generateTransparentSMatrix(matrixShape):
    STransparent = complexZeros((2, 2) + matrixShape);
    STransparent[0,1] = complexIdentity(matrixShape[0]);
    STransparent[1,0] = complexIdentity(matrixShape[0]);
    return STransparent;

def calculateRedhefferProduct(SA, SB):
    mat_shape = SA.shape;
    if(mat_shape != SB.shape):
        raise Exception(f'redhefferProduct: SA and SB are not of the same shape. SA is of shape {SA.shape} and SB is of shape {SB.shape}');

    SAB = complexZeros(mat_shape);
    D = calculateRedhefferDMatrix(SA, SB)
    F = calculateRedhefferFMatrix(SA, SB)

    SAB[0,0] = SA[0,0] + D @ SB[0,0] @ SA[1,0];
    SAB[0,1] = D @ SB[0,1];
    SAB[1,0] = F @ SA[1,0];
    SAB[1,1] = SB[1,1] +  F @ SA[1,1] @ SB[0,1];
    return SAB;

def calculatePMatrix(kx, ky, eri, uri):
    if isinstance(kx, np.ndarray):
        return calculatePMatrixNHarmonics(kx, ky, eri, uri)

    else:
        return calculatePMatrix1Harmonic(kx, ky, eri, uri)

def calculatePMatrix1Harmonic(kx, ky, eri, uri):
    P = complexZeros((2, 2));

    P[0,0] = kx*ky;
    P[0,1] = uri*eri - np.square(kx);
    P[1,0] = sq(ky) - uri*eri;
    P[1,1] = - kx*ky;
    P /= eri;
    return P

def calculatePMatrixNHarmonics(Kx, Ky, erConvolutionMatrix, urConvolutionMatrix):
    erConvolutionMatrixInverse = inv(erConvolutionMatrix)
    KMatrixDimension = Kx.shape[0]
    matrixShape = (2 *KMatrixDimension, 2 * KMatrixDimension)
    P = complexZeros(matrixShape)

    P[:KMatrixDimension,:KMatrixDimension] = Kx @ erConvolutionMatrixInverse @ Ky
    P[:KMatrixDimension,KMatrixDimension:] = urConvolutionMatrix - Kx @ erConvolutionMatrixInverse @ Kx
    P[KMatrixDimension:,:KMatrixDimension] = Ky @ erConvolutionMatrixInverse @ Ky - urConvolutionMatrix
    P[KMatrixDimension:,KMatrixDimension:] = - Ky @ erConvolutionMatrixInverse @ Kx
    return P

def calculateQMatrix(kx, ky, eri, uri):
    if isinstance(kx, np.ndarray):
        return calculateQMatrixNHarmonics(kx, ky, eri, uri)
    else:
        return calculateQMatrix1Harmonic(kx, ky, eri, uri)

def calculateQMatrix1Harmonic(kx, ky, eri, uri):
    Q = complexZeros(scatteringElementShape);

    Q[0,0] = kx * ky;
    Q[0,1] = uri*eri - sq(kx);
    Q[1,0] = sq(ky) - uri*eri;
    Q[1,1] = - kx * ky;
    Q = Q / uri;
    return Q;

def calculateQMatrixNHarmonics(Kx, Ky, erConvolutionMatrix, urConvolutionMatrix):
    urConvolutionMatrixInverse = inv(urConvolutionMatrix)
    KMatrixDimension = Kx.shape[0]
    matrixShape = (2 *KMatrixDimension, 2 * KMatrixDimension)
    Q = complexZeros(matrixShape)

    Q[:KMatrixDimension,:KMatrixDimension] = Kx @ urConvolutionMatrixInverse @ Ky
    Q[:KMatrixDimension,KMatrixDimension:] = erConvolutionMatrix - Kx @ urConvolutionMatrixInverse @ Kx
    Q[KMatrixDimension:,:KMatrixDimension] = Ky @ urConvolutionMatrixInverse @ Ky - erConvolutionMatrix
    Q[KMatrixDimension:,KMatrixDimension:] = - Ky @ urConvolutionMatrixInverse @ Kx
    return Q

def calculateOmegaMatrix(kz):
    return complexIdentity(2)* (0 + 1j)*kz;

def calculateOmegaSquaredMatrix(P, Q):
    return P @ Q

def calculateScatteringAMatrix(Wi, Wj, Vi, Vj):
    return inv(Wi) @ Wj + inv(Vi) @ Vj;

def calculateScatteringBMatrix(Wi, Wj, Vi, Vj): # UNIT TESTS COMPLETE
    return inv(Wi) @ Wj - inv(Vi) @ Vj;

def calculateScatteringDMatrix(Ai, Bi, Xi): # UNIT TESTS COMPLETE
    AiInverse = inv(Ai);
    return Ai - Xi @ Bi @ AiInverse @ Xi @ Bi;

def calculateRedhefferDMatrix(SA, SB):
    return SA[0,1] @ inv(complexIdentity(scatteringElementShape[0]) - SB[0,0] @ SA[1,1])

def calculateRedhefferFMatrix(SA, SB):
    return SB[1,0] @ inv(complexIdentity(scatteringElementShape[0]) - SB[0,0] @ SA[1,1])

def calculateKz(kx, ky, er, ur):
    return sqrt(er*ur - sq(kx) - sq(ky))

def calculateKzReflected(kx, ky, er, ur):
    return -calculateKz(kx, ky, er, ur)

def calculateKVector(theta, phi, er, ur):
    n = sqrt(er*ur);
    kx = n * sin(theta) * cos(phi);
    ky = n * sin(theta) * sin(phi);
    kz = n * cos(theta);
    return complexArray([kx, ky, kz]);

# This does not currently work because we don't have a Kz matrix in RCWA.
def calculateVWXMatrices(kx, ky, er, ur, k0=0, Li=0):
    if isinstance(kx, np.ndarray):
        return calculateVWXMatricesNHarmonics(kx, ky, er, ur, k0, Li)
    else:
        kz = calculateKz(kx, ky, er, ur)
        return calculateVWXMatrices1Harmonic(kx, ky, kz, er, ur, k0, Li)

def calculateVWXMatrices1Harmonic(kx, ky, kz, er, ur, k0, Li):
    Q = calculateQMatrix(kx, ky, er, ur);
    O = calculateOmegaMatrix(kz);
    OInverse = inv(O);
    W = complexIdentity(scatteringElementShape[0]);
    X = matrixExponentiate(O * k0 * Li)
    V = Q @ W @ OInverse;

    if(k0 > 0):
        return (V, W, X);
    else:
        return (V, W);

def calculateVWXMatricesNHarmonics(Kx, Ky, erConvolutionMatrix, urConvolutionMatrix, k0, Li):
    P = calculatePMatrix(Kx, Ky, erConvolutionMatrix, urConvolutionMatrix)
    Q = calculateQMatrix(Kx, Ky, erConvolutionMatrix, urConvolutionMatrix)
    OmegaSquared = calculateOmegaSquaredMatrix(P, Q)
    eigenValues, W = eig(OmegaSquared)
    Lambda = np.diag(np.sqrt(eigenValues))
    LambdaInverse = np.diag(np.reciprocal(sqrt(eigenValues)))
    V = Q @ W @ LambdaInverse
    X = matrixExponentiate(Lambda * k0 * Li)

    return (V, W, X)

def calculateInternalSMatrix(kx, ky, er, ur, k0, Li, Wg, Vg):
    (Vi, Wi, Xi) = calculateVWXMatrices(kx, ky, er, ur, k0, Li);
    Ai = calculateScatteringAMatrix(Wi, Wg, Vi, Vg);
    Bi = calculateScatteringBMatrix(Wi, Wg, Vi, Vg);
    Di = calculateScatteringDMatrix(Ai, Bi, Xi);

    Si = calculateInternalSMatrixFromRaw(Ai, Bi, Xi, Di);
    return Si;

def calculateReflectionRegionSMatrix(kx, ky, er, ur, Wg, Vg):
    kz = calculateKz(kx, ky, er, ur);
    (Vi, Wi) = calculateVWXMatrices(kx, ky, er, ur);
    Ai = calculateScatteringAMatrix(Wg, Wi, Vg, Vi);
    Bi = calculateScatteringBMatrix(Wg, Wi, Vg, Vi);

    Si = calculateReflectionRegionSMatrixFromRaw(Ai, Bi);
    return Si;

def calculateTransmissionRegionSMatrix(kx, ky, er, ur, Wg, Vg):
    (Vi, Wi) = calculateVWXMatrices(kx, ky, er, ur);
    Ai = calculateScatteringAMatrix(Wg, Wi, Vg, Vi);
    Bi = calculateScatteringBMatrix(Wg, Wi, Vg, Vi);

    Si = calculateTransmissionRegionSMatrixFromRaw(Ai, Bi);
    return Si;

def calculateInternalSMatrixFromRaw(Ai, Bi, Xi, Di):
    AiInverse = inv(Ai);
    DiInverse = inv(Di);

    S = complexZeros(scatteringMatrixShape);
    S[0,0] = DiInverse @ ((Xi @ Bi @ AiInverse @ Xi @ Ai) - Bi)
    S[0,1] = DiInverse @ Xi @ (Ai - (Bi @ AiInverse @ Bi));
    S[1,0] = S[0,1];
    S[1,1] = S[0,0];
    return S;

def calculateReflectionRegionSMatrixFromRaw(AReflectionRegion, BReflectionRegion):
    S = complexZeros(scatteringMatrixShape);
    A = AReflectionRegion;
    B = BReflectionRegion;
    AInverse = inv(A);

    S[0,0] = - AInverse @ B;
    S[0,1] = 2 * AInverse;
    S[1,0] = 0.5 * (A - B @ AInverse @ B)
    S[1,1] = B @ AInverse;
    return S;

def calculateTransmissionRegionSMatrixFromRaw(ATransmissionRegion, BTransmissionRegion): # UNIT TESTS COMPLETE
    A = ATransmissionRegion;
    B = BTransmissionRegion;
    AInverse = inv(A);

    S = complexZeros(scatteringMatrixShape);
    S[0,0] = B@ AInverse;
    S[0,1] = 0.5* (A- (B @ AInverse @ B))
    S[1,0] = 2* AInverse;
    S[1,1] = - AInverse @ B;
    return S;

def calculateEz(kx, ky, kz, Ex, Ey):
    Ez = - (kx*Ex + ky*Ey) / kz
    return Ez;

def calculateRT(kzReflectionRegion, kzTransmissionRegion,
        urReflectionRegion, urTransmissionRegion, ExyzReflected, ExyzTransmitted):
    R = sq(norm(ExyzReflected))
    T = sq(norm(ExyzTransmitted))*np.real(kzTransmissionRegion / urTransmissionRegion) / \
            (kzReflectionRegion / urReflectionRegion);

    return (R, T);
