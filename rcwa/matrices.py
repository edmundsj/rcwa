from rcwa.shorthand import *
from rcwa import Source, zeroSource
from numpy.linalg import pinv as pinv

def calculateIncidentFieldHarmonics(source, numberHarmonics):
    totalNumberHarmonics = np.prod(numberHarmonics)
    return np.hstack((source.pX * kroneckerDeltaVector(totalNumberHarmonics),
            source.pY * kroneckerDeltaVector(totalNumberHarmonics)))

def generateTransparentSMatrix(matrixShape):
    STransparent = complexZeros((2, 2) + matrixShape);
    STransparent[0,1] = complexIdentity(matrixShape[0]);
    STransparent[1,0] = complexIdentity(matrixShape[0]);
    return STransparent;

def calculateRedhefferProduct(SA, SB):
    SAB = complexZeros(SA.shape);
    D = calculateRedhefferDMatrix(SA, SB)
    F = calculateRedhefferFMatrix(SA, SB)

    SAB[0,0] = SA[0,0] + D @ SB[0,0] @ SA[1,0];
    SAB[0,1] = D @ SB[0,1];
    SAB[1,0] = F @ SA[1,0];
    SAB[1,1] = SB[1,1] +  F @ SA[1,1] @ SB[0,1];
    return SAB;

def calculatePMatrix(kx, ky, layer):
    if isinstance(kx, np.ndarray):
        return calculatePMatrixNHarmonics(kx, ky, layer)
    else:
        return calculatePMatrix1Harmonic(kx, ky, layer)

def calculatePMatrix1Harmonic(kx, ky, layer):
    P = complexZeros((2, 2));

    P[0,0] = kx*ky;
    P[0,1] = layer.er*layer.ur - np.square(kx);
    P[1,0] = sq(ky) - layer.er*layer.ur
    P[1,1] = - kx*ky;
    P /= layer.er;
    return P

def calculatePMatrixNHarmonics(Kx, Ky, layer):
    erInverse = pinv(layer.er)
    KMatrixDimension = Kx.shape[0]
    matrixShape = (2 *KMatrixDimension, 2 * KMatrixDimension)
    P = complexZeros(matrixShape)

    P[:KMatrixDimension,:KMatrixDimension] = Kx @ erInverse @ Ky
    P[:KMatrixDimension,KMatrixDimension:] = layer.ur - Kx @ erInverse @ Kx
    P[KMatrixDimension:,:KMatrixDimension] = Ky @ erInverse @ Ky - layer.ur
    P[KMatrixDimension:,KMatrixDimension:] = - Ky @ erInverse @ Kx
    return P

def calculateQMatrix(kx, ky, layer):
    if isinstance(kx, np.ndarray):
        if isinstance(layer.er, np.ndarray):
            return calculateQMatrixNHarmonics(kx, ky, layer)
        else:
            return calculateQReflectionTransmissionMatrix(kx, ky, layer)
    else:
        return calculateQMatrix1Harmonic(kx, ky, layer)

def calculateQMatrix1Harmonic(kx, ky, layer):
    Q = complexZeros((2,2));

    Q[0,0] = kx * ky;
    Q[0,1] = layer.er*layer.ur- sq(kx);
    Q[1,0] = sq(ky) - layer.er*layer.ur;
    Q[1,1] = - kx * ky;
    Q = Q / layer.ur;
    return Q;

def calculateQMatrixNHarmonics(Kx, Ky, layer):
    urInverse = pinv(layer.ur)
    KMatrixDimension = Kx.shape[0]
    matrixShape = (2 *KMatrixDimension, 2 * KMatrixDimension)
    Q = complexZeros(matrixShape)

    Q[:KMatrixDimension,:KMatrixDimension] = Kx @ urInverse @ Ky
    Q[:KMatrixDimension,KMatrixDimension:] = layer.er - Kx @ urInverse @ Kx
    Q[KMatrixDimension:,:KMatrixDimension] = Ky @ urInverse @ Ky - layer.er
    Q[KMatrixDimension:,KMatrixDimension:] = - Ky @ urInverse @ Kx
    return Q

def calculateQReflectionTransmissionMatrix(Kx, Ky, layer):
    KDimension = Kx.shape[0]
    Q = complexZeros((KDimension * 2, KDimension*2))
    Q[:KDimension, :KDimension] = Kx @ Ky
    Q[:KDimension, KDimension:] = layer.ur * layer.er *complexIdentity(KDimension) - Kx @ Kx
    Q[KDimension:, :KDimension] = Ky @ Ky - layer.ur*layer.er*complexIdentity(KDimension)
    Q[KDimension:, KDimension:] = - Ky @ Kx
    Q /= layer.ur
    return Q

def calculateLambdaMatrix(kz):
    if isinstance(kz, np.ndarray):
        KzDimension = kz.shape[0]
        LambdaShape = (KzDimension*2, KzDimension*2)
        Lambda = complexZeros(LambdaShape)
        Lambda[:KzDimension, :KzDimension] = 1j*kz
        Lambda[KzDimension:, KzDimension:] = 1j*kz
        return Lambda
    else:
        return complexIdentity(2)* (0 + 1j)*kz;

def calculateOmegaSquaredMatrix(P, Q):
    return P @ Q

def calculateScatteringAMatrix(Wi, Wj, Vi, Vj):
    return pinv(Wi) @ Wj + inv(Vi) @ Vj;

def calculateScatteringBMatrix(Wi, Wj, Vi, Vj):
    return pinv(Wi) @ Wj - inv(Vi) @ Vj;

def calculateScatteringDMatrix(Ai, Bi, Xi):
    AiInverse = pinv(Ai);
    return Ai - Xi @ Bi @ AiInverse @ Xi @ Bi;

def calculateRedhefferDMatrix(SA, SB):
    return SA[0,1] @ pinv(complexIdentity(SA[0,0].shape[0]) - SB[0,0] @ SA[1,1])

def calculateRedhefferFMatrix(SA, SB):
    return SB[1,0] @ pinv(complexIdentity(SA[0,0].shape[0]) - SA[1,1] @ SB[0,0])

def calculateKzBackward(kx, ky, layer):
    if isinstance(kx, np.ndarray):
        return -conj(sqrt(conj(layer.er*layer.ur)*complexIdentity(kx.shape[0]) - kx @ kx - ky @ ky))
    else:
        return sqrt(layer.er*layer.ur - sq(kx) - sq(ky))

def calculateKzForward(kx, ky, layer):
    if isinstance(kx, np.ndarray):
        return conj(sqrt(conj(layer.er*layer.ur)*complexIdentity(kx.shape[0]) - kx @ kx - ky @ ky))
    else:
        return sqrt(layer.er*layer.ur - sq(kx) - sq(ky))

def calculateKVector(source, layer):
    kx = layer.n * sin(source.theta) * cos(source.phi);
    ky = layer.n * sin(source.theta) * sin(source.phi);
    kz = layer.n * cos(source.theta);
    return complexArray([kx, ky, kz]);

# TODO: MUST BE MODIFIED TO HANDLE HOMOGENOUS LAYERS AND NOT DIRECTLY SOLVE THE 
# EIGENVALUE PROBLEM, BUT PROVIDE THE KNOWN SOLUTION.
def calculateVWXMatrices(kx, ky, layer, source=zeroSource):
    if isinstance(kx, np.ndarray):
        return calculateVWXMatricesNHarmonics(kx, ky, layer, source)
    else:
        kz = calculateKzForward(kx, ky, layer) # This *should* work.
        return calculateVWXMatrices1Harmonic(kx, ky, kz, layer, source)

def calculateVWXMatrices1Harmonic(kx, ky, kz, layer, source):
    Q = calculateQMatrix(kx, ky, layer);
    O = calculateLambdaMatrix(kz);
    OInverse = pinv(O);
    W = complexIdentity(2)
    X = matrixExponentiate(O * source.k0 * layer.L)
    V = Q @ W @ OInverse;

    return (V, W, X);

# TODO - INCORPORATE CONVOLUTION MATRIX INTO THE LAYERS
def calculateVWXMatricesNHarmonics(Kx, Ky, layer, source):
    P = calculatePMatrix(Kx, Ky, layer)
    Q = calculateQMatrix(Kx, Ky, layer)
    OmegaSquared = calculateOmegaSquaredMatrix(P, Q)

    if layer.homogenous is False:
        eigenValues, W = eig(OmegaSquared)
        Lambda = np.diag(sqrt(eigenValues))
        LambdaInverse = np.diag(np.reciprocal(sqrt(eigenValues)))
        V = Q @ W @ LambdaInverse
        X = matrixExponentiate( -Lambda * source.k0 * layer.L)
        return (V, W, X)
    else:
        Kz = calculateKzForward(Kx, Ky, layer)
        Lambda = calculateLambdaMatrix(Kz)
        LambdaInverse = pinv(Lambda)
        W = complexIdentity(2*Kz.shape[0])
        V = Q @ W @ LambdaInverse
        X = matrixExponentiate( -Lambda * source.k0 * layer.L)
        return (V, W, X)

def calculateInternalSMatrix(kx, ky, layer, source, Wg, Vg):
    (Vi, Wi, Xi) = calculateVWXMatrices(kx, ky, layer, source);
    Ai = calculateScatteringAMatrix(Wi, Wg, Vi, Vg);
    Bi = calculateScatteringBMatrix(Wi, Wg, Vi, Vg);
    Di = calculateScatteringDMatrix(Ai, Bi, Xi);

    Si = calculateInternalSMatrixFromRaw(Ai, Bi, Xi, Di);
    return Si;

def calculateReflectionRegionSMatrix(kx, ky, layerStack, Wg, Vg):
    if isinstance(kx, np.ndarray):
        return calculateReflectionRegionSMatrixNHarmonics(kx, ky, layerStack, Wg, Vg)
    else:
        return calculateReflectionRegionSMatrix1Harmonic(kx, ky, layerStack, Wg, Vg)

def calculateReflectionRegionSMatrix1Harmonic(kx, ky, layerStack, Wg, Vg):
    kz = calculateKzForward(kx, ky, layerStack.reflectionLayer);
    (Vi, Wi, X) = calculateVWXMatrices(kx, ky, layerStack.reflectionLayer);
    Ai = calculateScatteringAMatrix(Wg, Wi, Vg, Vi);
    Bi = calculateScatteringBMatrix(Wg, Wi, Vg, Vi);

    Si = calculateReflectionRegionSMatrixFromRaw(Ai, Bi);
    return Si;

def calculateReflectionRegionSMatrixNHarmonics(Kx, Ky, layerStack, Wg, Vg):
    KDimension = Kx.shape[0]
    lambdaRef = complexZeros((KDimension*2, KDimension*2))
    Wi = complexIdentity(KDimension * 2)
    Q = calculateQMatrix(Kx, Ky, layerStack.reflectionLayer)

    # I have no idea why we conjugate ur * er and then conjugate the whole thing.
    Kz = conj(sqrt (conj(layerStack.reflectionLayer.er*layerStack.reflectionLayer.ur)* \
            complexIdentity(KDimension) - Kx @ Kx - Ky @ Ky))
    lambdaRef[:KDimension, :KDimension] = 1j*Kz
    lambdaRef[KDimension:,KDimension:] = 1j*Kz
    Vi = Q @ pinv(lambdaRef)
    Ai = calculateScatteringAMatrix(Wg, Wi, Vg, Vi)
    Bi = calculateScatteringBMatrix(Wg, Wi, Vg, Vi)

    Sref = calculateReflectionRegionSMatrixFromRaw(Ai, Bi)
    return Sref

def calculateTransmissionRegionSMatrix(kx, ky, layerStack, Wg, Vg):
    if isinstance(kx, np.ndarray):
        return calculateTransmissionRegionSMatrixNHarmonics(kx, ky, layerStack, Wg, Vg)
    else:
        return calculateTransmissionRegionSMatrix1Harmonic(kx, ky, layerStack, Wg, Vg)

def calculateTransmissionRegionSMatrix1Harmonic(kx, ky, layerStack, Wg, Vg):
    (Vi, Wi, X) = calculateVWXMatrices(kx, ky, layerStack.transmissionLayer)
    Ai = calculateScatteringAMatrix(Wg, Wi, Vg, Vi);
    Bi = calculateScatteringBMatrix(Wg, Wi, Vg, Vi);

    Si = calculateTransmissionRegionSMatrixFromRaw(Ai, Bi);
    return Si;

def calculateTransmissionRegionSMatrixNHarmonics(Kx, Ky, layerStack, Wg, Vg):
    KDimension = Kx.shape[0]
    lambdaRef = complexZeros((KDimension*2, KDimension*2))
    Wi = complexIdentity(KDimension * 2)
    Q = calculateQMatrix(Kx, Ky, layerStack.transmissionLayer)

    # I have no idea why we conjugate ur * er and then conjugate the whole thing.
    Kz = conj(sqrt (conj(layerStack.transmissionLayer.er * layerStack.transmissionLayer.ur)* \
            complexIdentity(KDimension) - Kx @ Kx - Ky @ Ky))
    lambdaRef[:KDimension, :KDimension] = 1j*Kz
    lambdaRef[KDimension:,KDimension:] = 1j*Kz
    Vi = Q @ pinv(lambdaRef)
    Ai = calculateScatteringAMatrix(Wg, Wi, Vg, Vi)
    Bi = calculateScatteringBMatrix(Wg, Wi, Vg, Vi)

    Strn = calculateTransmissionRegionSMatrixFromRaw(Ai, Bi)
    return Strn

def calculateInternalSMatrixFromRaw(Ai, Bi, Xi, Di):
    AiInverse = pinv(Ai)
    DiInverse = pinv(Di);

    S = complexZeros((2, 2) + Ai.shape)
    S[0, 0] = DiInverse @ (Xi @ Bi @ AiInverse @ Xi @ Ai - Bi)
    S[0, 1] = DiInverse @ Xi @ (Ai - Bi @ AiInverse @ Bi)
    S[1, 0] = S[0, 1]
    S[1, 1] = S[0, 0]
    return S

def calculateReflectionRegionSMatrixFromRaw(AReflectionRegion, BReflectionRegion):
    S = complexZeros((2, 2) + AReflectionRegion.shape);
    A = AReflectionRegion;
    B = BReflectionRegion;
    AInverse = pinv(A);

    S[0,0] = - AInverse @ B;
    S[0,1] = 2 * AInverse;
    S[1,0] = 0.5 * (A - B @ AInverse @ B)
    S[1,1] = B @ AInverse;
    return S;

def calculateTransmissionRegionSMatrixFromRaw(ATransmissionRegion, BTransmissionRegion): # UNIT TESTS COMPLETE
    A = ATransmissionRegion;
    B = BTransmissionRegion;
    AInverse = pinv(A);

    S = complexZeros((2, 2) + A.shape);
    S[0,0] = B@ AInverse;
    S[0,1] = 0.5* (A- (B @ AInverse @ B))
    S[1,0] = 2* AInverse;
    S[1,1] = - AInverse @ B;
    return S;

# NOTE - this can only be used for 1D (TMM-type) simulations. rTE/rTM are not meaningful quantities otherwise.
def calculateTEMReflectionCoefficientsFromXYZ(source, rx, ry, rz):
    if isinstance(rx, np.ndarray):
        raise NotImplementedError
    else:
        rxyz = np.array([rx, ry, rz])
        rTEM = source.ATEM @ rxyz
        rTEM[0] = rTEM[0] / source.pTE
        rTEM[1] = rTEM[1] / source.pTM
        return rTEM

def calculateReflectionCoefficient(S, Kx, Ky, KzReflectionRegion,
        WReflectionRegion, source, numberHarmonics):
    incidentFieldHarmonics = calculateIncidentFieldHarmonics(source, numberHarmonics)
    rTransverse = WReflectionRegion @ S[0,0] @ pinv(WReflectionRegion) @ incidentFieldHarmonics

    rx, ry, rz = None, None, None
    if isinstance(Kx, np.ndarray):
        maxIndex = int(rTransverse.shape[0]/2)
        rx = rTransverse[0:maxIndex]
        ry = rTransverse[maxIndex:]
        rz = - inv(KzReflectionRegion) @ (Kx @ rx + Ky @ ry)
    else:
        rx = rTransverse[0]
        ry = rTransverse[1]
        rz = - (Kx * rx + Ky * ry) / KzReflectionRegion
    return rx, ry, rz

def calculateTransmissionCoefficient(S, Kx, Ky, KzTransmissionRegion,
        WTransmissionRegion, source, numberHarmonics):
    incidentFieldHarmonics = calculateIncidentFieldHarmonics(source, numberHarmonics)
    tTransverse = WTransmissionRegion @ S[1,0] @ inv(WTransmissionRegion) @ incidentFieldHarmonics

    tx, ty, tz = None, None, None
    if isinstance(Kx, np.ndarray):
        maxIndex = int(tTransverse.shape[0]/2)
        tx = tTransverse[:maxIndex]
        ty = tTransverse[maxIndex:]
        tz = - inv(KzTransmissionRegion) @ (Kx @ tx + Ky @ ty)
    else:
        tx = tTransverse[0]
        ty = tTransverse[1]
        tz = - (Kx * tx + Ky * ty) / KzTransmissionRegion
    return tx, ty, tz

def calculateDiffractionReflectionEfficiency(rx, ry, rz, source, KzReflectionRegion, layerStack):
    urReflectionRegion = layerStack.reflectionLayer.ur
    preMatrix = real(-1 /urReflectionRegion * KzReflectionRegion) / \
            real(source.kIncident[2] / urReflectionRegion)
    R = None
    if isinstance(KzReflectionRegion, np.ndarray):
        R = preMatrix @ (sq(np.abs(rx)) + sq(np.abs(ry)) + sq(np.abs(rz)))
        RDimension = int(sqrt(rx.shape[0]))
        R = R.reshape((RDimension, RDimension))
    else:
        R = -preMatrix * (sq(np.abs(rx)) + sq(np.abs(ry)) + sq(np.abs(rz)))
    return R

def calculateDiffractionTransmissionEfficiency(tx, ty, tz, source, KzTransmissionRegion, layerStack):
    urTransmissionRegion = layerStack.transmissionLayer.ur
    urReflectionRegion = layerStack.reflectionLayer.ur
    preMatrix = real(1 / urTransmissionRegion * KzTransmissionRegion) / \
            real(source.kIncident[2] / urReflectionRegion)

    if isinstance(KzTransmissionRegion, np.ndarray):
        T = preMatrix @ (sq(np.abs(tx)) + sq(np.abs(ty)) + sq(np.abs(tz)))
        TDimension = int(sqrt(tx.shape[0]))
        T = T.reshape((TDimension, TDimension))
    else:
        T = preMatrix * (sq(np.abs(tx)) + sq(np.abs(ty)) + sq(np.abs(tz)))
    return T

def calculateEz(kx, ky, kz, Ex, Ey):
    Ez = - (kx*Ex + ky*Ey) / kz
    return Ez;

def calculateRT(kzReflectionRegion, kzTransmissionRegion,
        layerStack, ExyzReflected, ExyzTransmitted):
    urTransmissionRegion = layerStack.transmissionLayer.ur
    urReflectionRegion = layerStack.reflectionLayer.ur
    R = sq(norm(ExyzReflected))
    T = sq(norm(ExyzTransmitted))*np.real(kzTransmissionRegion / urTransmissionRegion) / \
            (kzReflectionRegion / urReflectionRegion);

    return (R, T);
