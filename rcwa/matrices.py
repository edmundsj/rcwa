from rcwa.shorthand import *
from rcwa import Source, zeroSource
from numpy.linalg import pinv as pinv

def s_incident(source, numberHarmonics):
    totalNumberHarmonics = np.prod(numberHarmonics)
    return np.hstack((source.pX * kroneckerDeltaVector(totalNumberHarmonics),
            source.pY * kroneckerDeltaVector(totalNumberHarmonics)))

def S_matrix_transparent(matrixShape):
    STransparent = complexZeros((2, 2) + matrixShape);
    STransparent[0,1] = complexIdentity(matrixShape[0]);
    STransparent[1,0] = complexIdentity(matrixShape[0]);
    return STransparent;

def redheffer_product(SA, SB):
    SAB = complexZeros(SA.shape);
    D = D_matrix_redheffer(SA, SB)
    F = F_matrix(SA, SB)

    SAB[0,0] = SA[0,0] + D @ SB[0,0] @ SA[1,0];
    SAB[0,1] = D @ SB[0,1];
    SAB[1,0] = F @ SA[1,0];
    SAB[1,1] = SB[1,1] +  F @ SA[1,1] @ SB[0,1];
    return SAB;

def P_matrix(kx, ky, layer):
    if isinstance(kx, np.ndarray):
        return P_matrix_general(kx, ky, layer)
    else:
        return P_matrix_homogenous(kx, ky, layer)

def P_matrix_homogenous(kx, ky, layer):
    P = complexZeros((2, 2));

    P[0,0] = kx*ky;
    P[0,1] = layer.er*layer.ur - np.square(kx);
    P[1,0] = sq(ky) - layer.er*layer.ur
    P[1,1] = - kx*ky;
    P /= layer.er;
    return P

def P_matrix_general(Kx, Ky, layer):
    erInverse = pinv(layer.er)
    KMatrixDimension = Kx.shape[0]
    matrixShape = (2 *KMatrixDimension, 2 * KMatrixDimension)
    P = complexZeros(matrixShape)

    P[:KMatrixDimension,:KMatrixDimension] = Kx @ erInverse @ Ky
    P[:KMatrixDimension,KMatrixDimension:] = layer.ur - Kx @ erInverse @ Kx
    P[KMatrixDimension:,:KMatrixDimension] = Ky @ erInverse @ Ky - layer.ur
    P[KMatrixDimension:,KMatrixDimension:] = - Ky @ erInverse @ Kx
    return P

def Q_matrix(kx, ky, layer):
    if isinstance(kx, np.ndarray):
        if isinstance(layer.er, np.ndarray):
            return Q_matrix_general(kx, ky, layer)
        else:
            return Q_matrix_semi_infinite(kx, ky, layer)
    else:
        return Q_matrix_homogenous(kx, ky, layer)

def Q_matrix_homogenous(kx, ky, layer):
    Q = complexZeros((2,2));

    Q[0,0] = kx * ky;
    Q[0,1] = layer.er*layer.ur- sq(kx);
    Q[1,0] = sq(ky) - layer.er*layer.ur;
    Q[1,1] = - kx * ky;
    Q = Q / layer.ur;
    return Q;

def Q_matrix_general(Kx, Ky, layer):
    urInverse = pinv(layer.ur)
    KMatrixDimension = Kx.shape[0]
    matrixShape = (2 *KMatrixDimension, 2 * KMatrixDimension)
    Q = complexZeros(matrixShape)

    Q[:KMatrixDimension,:KMatrixDimension] = Kx @ urInverse @ Ky
    Q[:KMatrixDimension,KMatrixDimension:] = layer.er - Kx @ urInverse @ Kx
    Q[KMatrixDimension:,:KMatrixDimension] = Ky @ urInverse @ Ky - layer.er
    Q[KMatrixDimension:,KMatrixDimension:] = - Ky @ urInverse @ Kx
    return Q

def Q_matrix_semi_infinite(Kx, Ky, layer):
    KDimension = Kx.shape[0]
    Q = complexZeros((KDimension * 2, KDimension*2))
    Q[:KDimension, :KDimension] = Kx @ Ky
    Q[:KDimension, KDimension:] = layer.ur * layer.er *complexIdentity(KDimension) - Kx @ Kx
    Q[KDimension:, :KDimension] = Ky @ Ky - layer.ur*layer.er*complexIdentity(KDimension)
    Q[KDimension:, KDimension:] = - Ky @ Kx
    Q /= layer.ur
    return Q

def lambda_matrix(kz):
    if isinstance(kz, np.ndarray):
        KzDimension = kz.shape[0]
        LambdaShape = (KzDimension*2, KzDimension*2)
        Lambda = complexZeros(LambdaShape)
        Lambda[:KzDimension, :KzDimension] = 1j*kz
        Lambda[KzDimension:, KzDimension:] = 1j*kz
        return Lambda
    else:
        return complexIdentity(2)* (0 + 1j)*kz;

def omega_squared_matrix(P, Q):
    return P @ Q

def A_matrix(Wi, Wj, Vi, Vj):
    return pinv(Wi) @ Wj + inv(Vi) @ Vj;

def B_matrix(Wi, Wj, Vi, Vj):
    return pinv(Wi) @ Wj - inv(Vi) @ Vj;

def D_matrix(Ai, Bi, Xi):
    AiInverse = pinv(Ai);
    return Ai - Xi @ Bi @ AiInverse @ Xi @ Bi;

def D_matrix_redheffer(SA, SB):
    return SA[0,1] @ pinv(complexIdentity(SA[0,0].shape[0]) - SB[0,0] @ SA[1,1])

def F_matrix(SA, SB):
    return SB[1,0] @ pinv(complexIdentity(SA[0,0].shape[0]) - SA[1,1] @ SB[0,0])

def Kz_backward(kx, ky, layer):
    if isinstance(kx, np.ndarray):
        return -conj(sqrt(conj(layer.er*layer.ur)*complexIdentity(kx.shape[0]) - kx @ kx - ky @ ky))
    else:
        return sqrt(layer.er*layer.ur - sq(kx) - sq(ky))

def Kz_forward(kx, ky, layer):
    if isinstance(kx, np.ndarray):
        return conj(sqrt(conj(layer.er*layer.ur)*complexIdentity(kx.shape[0]) - kx @ kx - ky @ ky))
    else:
        return sqrt(layer.er*layer.ur - sq(kx) - sq(ky))


def k_vector(source, layer):
    kx = layer.n * sin(source.theta) * cos(source.phi);
    ky = layer.n * sin(source.theta) * sin(source.phi);
    kz = layer.n * cos(source.theta);
    return complexArray([kx, ky, kz]);


def VWX_matrices(kx, ky, layer, source=zeroSource):
    if not isinstance(kx, np.ndarray):
        kz = Kz_forward(kx, ky, layer)
        return VWX_matrices_homogenous(kx, ky, kz, layer, source)
    else:
        return VWX_matrices_general(kx, ky, layer, source)


def VWX_matrices_homogenous(kx, ky, kz, layer, source):
    Q = Q_matrix(kx, ky, layer);
    O = lambda_matrix(kz);
    OInverse = pinv(O);
    W = complexIdentity(2)
    X = matrixExponentiate(O * source.k0 * layer.thickness)
    V = Q @ W @ OInverse;

    return (V, W, X);


def VWX_matrices_general(Kx, Ky, layer, source):
    P = P_matrix(Kx, Ky, layer)
    Q = Q_matrix(Kx, Ky, layer)
    OmegaSquared = omega_squared_matrix(P, Q)

    if layer.homogenous:
        Kz = Kz_forward(Kx, Ky, layer)
        Lambda = lambda_matrix(Kz)
        LambdaInverse = pinv(Lambda)
        W = complexIdentity(2 * Kz.shape[0])
        V = Q @ W @ LambdaInverse
        X = matrixExponentiate(-Lambda * source.k0 * layer.thickness)
        return (V, W, X)
    else:
        eigenValues, W = eig(OmegaSquared)
        Lambda = np.diag(sqrt(eigenValues))
        LambdaInverse = np.diag(np.reciprocal(sqrt(eigenValues)))
        V = Q @ W @ LambdaInverse
        X = matrixExponentiate(-Lambda * source.k0 * layer.thickness)
        return (V, W, X)

def S_matrix_internal(kx, ky, layer, source, Wg, Vg):
    (Vi, Wi, Xi) = VWX_matrices(kx, ky, layer, source);
    Ai = A_matrix(Wi, Wg, Vi, Vg);
    Bi = B_matrix(Wi, Wg, Vi, Vg);
    Di = D_matrix(Ai, Bi, Xi);

    Si = calculateInternalSMatrixFromRaw(Ai, Bi, Xi, Di);
    return Si;

def calculateReflectionRegionSMatrix(kx, ky, layerStack, Wg, Vg):
    if isinstance(kx, np.ndarray):
        return calculateReflectionRegionSMatrixNHarmonics(kx, ky, layerStack, Wg, Vg)
    else:
        return calculateReflectionRegionSMatrix1Harmonic(kx, ky, layerStack, Wg, Vg)

def calculateReflectionRegionSMatrix1Harmonic(kx, ky, layerStack, Wg, Vg):
    kz = Kz_forward(kx, ky, layerStack.incident_layer);
    (Vi, Wi, X) = VWX_matrices(kx, ky, layerStack.incident_layer);
    Ai = A_matrix(Wg, Wi, Vg, Vi);
    Bi = B_matrix(Wg, Wi, Vg, Vi);

    Si = calculateReflectionRegionSMatrixFromRaw(Ai, Bi);
    return Si;

def calculateReflectionRegionSMatrixNHarmonics(Kx, Ky, layerStack, Wg, Vg):
    KDimension = Kx.shape[0]
    lambdaRef = complexZeros((KDimension*2, KDimension*2))
    Wi = complexIdentity(KDimension * 2)
    Q = Q_matrix(Kx, Ky, layerStack.incident_layer)

    # I have no idea why we conjugate ur * er and then conjugate the whole thing.
    Kz = conj(sqrt (conj(layerStack.incident_layer.er * layerStack.incident_layer.ur) * \
                    complexIdentity(KDimension) - Kx @ Kx - Ky @ Ky))
    lambdaRef[:KDimension, :KDimension] = 1j*Kz
    lambdaRef[KDimension:,KDimension:] = 1j*Kz
    Vi = Q @ pinv(lambdaRef)
    Ai = A_matrix(Wg, Wi, Vg, Vi)
    Bi = B_matrix(Wg, Wi, Vg, Vi)

    Sref = calculateReflectionRegionSMatrixFromRaw(Ai, Bi)
    return Sref

def calculateTransmissionRegionSMatrix(kx, ky, layerStack, Wg, Vg):
    if isinstance(kx, np.ndarray):
        return calculateTransmissionRegionSMatrixNHarmonics(kx, ky, layerStack, Wg, Vg)
    else:
        return calculateTransmissionRegionSMatrix1Harmonic(kx, ky, layerStack, Wg, Vg)

def calculateTransmissionRegionSMatrix1Harmonic(kx, ky, layerStack, Wg, Vg):
    (Vi, Wi, X) = VWX_matrices(kx, ky, layerStack.transmission_layer)
    Ai = A_matrix(Wg, Wi, Vg, Vi);
    Bi = B_matrix(Wg, Wi, Vg, Vi);

    Si = calculateTransmissionRegionSMatrixFromRaw(Ai, Bi);
    return Si;

def calculateTransmissionRegionSMatrixNHarmonics(Kx, Ky, layerStack, Wg, Vg):
    KDimension = Kx.shape[0]
    lambdaRef = complexZeros((KDimension*2, KDimension*2))
    Wi = complexIdentity(KDimension * 2)
    Q = Q_matrix(Kx, Ky, layerStack.transmission_layer)

    # I have no idea why we conjugate ur * er and then conjugate the whole thing.
    Kz = conj(sqrt (conj(layerStack.transmission_layer.er * layerStack.transmission_layer.ur) * \
                    complexIdentity(KDimension) - Kx @ Kx - Ky @ Ky))
    lambdaRef[:KDimension, :KDimension] = 1j*Kz
    lambdaRef[KDimension:,KDimension:] = 1j*Kz
    Vi = Q @ pinv(lambdaRef)
    Ai = A_matrix(Wg, Wi, Vg, Vi)
    Bi = B_matrix(Wg, Wi, Vg, Vi)

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
    incidentFieldHarmonics = s_incident(source, numberHarmonics)
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
    incidentFieldHarmonics = s_incident(source, numberHarmonics)
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

def calculateDiffractionReflectionEfficiency(rx, ry, rz, source, KzReflectionRegion, layerStack, numberHarmonics):
    urReflectionRegion = layerStack.incident_layer.ur
    preMatrix = real(-1 /urReflectionRegion * KzReflectionRegion) / \
            real(source.k_incident[2] / urReflectionRegion)
    R = None
    if isinstance(KzReflectionRegion, np.ndarray):
        R = preMatrix @ (sq(np.abs(rx)) + sq(np.abs(ry)) + sq(np.abs(rz)))
        RDimension = int(sqrt(rx.shape[0]))
        if not np.isscalar(numberHarmonics):
            R = R.reshape((RDimension, RDimension))
    else:
        R = -preMatrix * (sq(np.abs(rx)) + sq(np.abs(ry)) + sq(np.abs(rz)))
    return R

def calculateDiffractionTransmissionEfficiency(tx, ty, tz, source, KzTransmissionRegion, layerStack,
                                              numberHarmonics):
    urTransmissionRegion = layerStack.transmission_layer.ur
    urReflectionRegion = layerStack.incident_layer.ur
    preMatrix = real(1 / urTransmissionRegion * KzTransmissionRegion) / \
            real(source.k_incident[2] / urReflectionRegion)

    if isinstance(KzTransmissionRegion, np.ndarray):
        T = preMatrix @ (sq(np.abs(tx)) + sq(np.abs(ty)) + sq(np.abs(tz)))
        TDimension = int(sqrt(tx.shape[0]))
        if not np.isscalar(numberHarmonics):
            T = T.reshape((TDimension, TDimension))
    else:
        T = preMatrix * (sq(np.abs(tx)) + sq(np.abs(ty)) + sq(np.abs(tz)))
    return T

def calculateEz(kx, ky, kz, Ex, Ey):
    Ez = - (kx*Ex + ky*Ey) / kz
    return Ez;

def calculateRT(kzReflectionRegion, kzTransmissionRegion,
        layerStack, ExyzReflected, ExyzTransmitted):
    urTransmissionRegion = layerStack.transmission_layer.ur
    urReflectionRegion = layerStack.incident_layer.ur
    R = sq(norm(ExyzReflected))
    T = sq(norm(ExyzTransmitted))*np.real(kzTransmissionRegion / urTransmissionRegion) / \
            (kzReflectionRegion / urReflectionRegion);

    return (R, T);
