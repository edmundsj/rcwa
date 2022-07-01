from rcwa.shorthand import *
from numpy.linalg import pinv as pinv

def s_incident(source, n_harmonics):
    totalNumberHarmonics = np.prod(n_harmonics)
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

def k_vector(source, layer):
    kx = layer.n * sin(source.theta) * cos(source.phi)
    ky = layer.n * sin(source.theta) * sin(source.phi)
    kz = layer.n * cos(source.theta)
    return complexArray([kx, ky, kz])

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


class MatrixCalculator:
    """
    Superclass of Layer which is used purely for the calculation of matrices
    """

    def P_matrix(self):
        if isinstance(self.Kx, np.ndarray):
            return self._P_matrix_general()
        else:
            return self._P_matrix_homogenous()

    def _P_matrix_homogenous(self):
        P = complexZeros((2, 2));

        P[0,0] = self.Kx*self.Ky;
        P[0,1] = self.er*self.ur - np.square(self.Kx);
        P[1,0] = sq(self.Ky) - self.er*self.ur
        P[1,1] = - self.Kx*self.Ky;
        P /= self.er;
        return P

    def _P_matrix_general(self):
        erInverse = pinv(self.er)
        KMatrixDimension = self.Kx.shape[0]
        matrixShape = (2 *KMatrixDimension, 2 * KMatrixDimension)
        P = complexZeros(matrixShape)

        P[:KMatrixDimension,:KMatrixDimension] = self.Kx @ erInverse @ self.Ky
        P[:KMatrixDimension,KMatrixDimension:] = self.ur - self.Kx @ erInverse @ self.Kx
        P[KMatrixDimension:,:KMatrixDimension] = self.Ky @ erInverse @ self.Ky - self.ur
        P[KMatrixDimension:,KMatrixDimension:] = - self.Ky @ erInverse @ self.Kx
        return P

    def Q_matrix(self):
        if isinstance(self.Kx, np.ndarray):
            if isinstance(self.er, np.ndarray):
                return self._Q_matrix_general()
            else:
                return self._Q_matrix_semi_infinite()
        else:
            return self._Q_matrix_homogenous()

    def _Q_matrix_homogenous(self):
        Q = complexZeros((2,2));

        Q[0,0] = self.Kx * self.Ky;
        Q[0,1] = self.er*self.ur - sq(self.Kx);
        Q[1,0] = sq(self.Ky) - self.er*self.ur;
        Q[1,1] = - self.Kx * self.Ky;
        Q = Q / self.ur;
        return Q;

    def _Q_matrix_general(self):
        urInverse = pinv(self.ur)
        KMatrixDimension = self.Kx.shape[0]
        matrixShape = (2 *KMatrixDimension, 2 * KMatrixDimension)
        Q = complexZeros(matrixShape)

        Q[:KMatrixDimension,:KMatrixDimension] = self.Kx @ urInverse @ self.Ky
        Q[:KMatrixDimension,KMatrixDimension:] = self.er - self.Kx @ urInverse @ self.Kx
        Q[KMatrixDimension:,:KMatrixDimension] = self.Ky @ urInverse @ self.Ky - self.er
        Q[KMatrixDimension:,KMatrixDimension:] = - self.Ky @ urInverse @ self.Kx
        return Q

    def _Q_matrix_semi_infinite(self):
        KDimension = self.Kx.shape[0]
        Q = complexZeros((KDimension * 2, KDimension*2))
        Q[:KDimension, :KDimension] = self.Kx @ self.Ky
        Q[:KDimension, KDimension:] = self.ur * self.er * complexIdentity(KDimension) - self.Kx @ self.Kx
        Q[KDimension:, :KDimension] = self.Ky @ self.Ky - self.ur*self.er*complexIdentity(KDimension)
        Q[KDimension:, KDimension:] = - self.Ky @ self.Kx
        Q /= self.ur
        return Q

    def lambda_matrix(self):
        Kz = self.Kz_forward() # I am a little unsure about this particular line. Why is Kz_backward never used?

        if isinstance(Kz, np.ndarray):
            KzDimension = Kz.shape[0]
            LambdaShape = (KzDimension*2, KzDimension*2)
            Lambda = complexZeros(LambdaShape)
            Lambda[:KzDimension, :KzDimension] = 1j*Kz
            Lambda[KzDimension:, KzDimension:] = 1j*Kz
            return Lambda
        else:
            return complexIdentity(2)* (0 + 1j)*Kz;

    def Kz_backward(self):
        if isinstance(self.Kx, np.ndarray):
            return -conj(sqrt(conj(self.er*self.ur)*complexIdentity(self.Kx.shape[0]) - self.Kx @ self.Kx - self.Ky @ self.Ky))
        else:
            return sqrt(self.er*self.ur - sq(self.Kx) - sq(self.Ky))

    def Kz_forward(self):
        if isinstance(self.Kx, np.ndarray):
            return conj(sqrt(conj(self.er*self.ur)*complexIdentity(self.Kx.shape[0]) - self.Kx @ self.Kx - self.Ky @ self.Ky))
        else:
            return sqrt(self.er*self.ur - sq(self.Kx) - sq(self.Ky))

    def Kz_gap(self):
        if isinstance(self.Kx, np.ndarray):
            return conj(sqrt(complexIdentity(self.Kx.shape[0]) - self.Kx @ self.Kx - self.Ky @ self.Ky))
        else:
            return sqrt(self.er*self.ur - sq(self.Kx) - sq(self.Ky))

    def VWX_matrices(self):
        if not isinstance(self.Kx, np.ndarray):
            return self._VWX_matrices_homogenous()
        else:
            return self._VWX_matrices_general()

    def _VWX_matrices_homogenous(self):
        Kz = self.Kz_forward()
        Q = self.Q_matrix()
        O = self.lambda_matrix()
        OInverse = pinv(O)
        W = complexIdentity(2)
        X = matrixExponentiate(O * self.source.k0 * self.thickness)
        V = Q @ W @ OInverse

        return (V, W, X)

    def _VWX_matrices_general(self):
        P = self.P_matrix()
        Q = self.Q_matrix()
        OmegaSquared = omega_squared_matrix(P, Q)

        if self.homogenous:
            Kz = self.Kz_forward()
            Lambda = self.lambda_matrix()
            LambdaInverse = pinv(Lambda)
            W = complexIdentity(2 * Kz.shape[0])
            V = Q @ W @ LambdaInverse
            X = matrixExponentiate(-Lambda * self.source.k0 * self.thickness)
            return (V, W, X)
        else:
            eigenValues, W = eig(OmegaSquared)
            Lambda = np.diag(sqrt(eigenValues))
            LambdaInverse = np.diag(np.reciprocal(sqrt(eigenValues)))
            V = Q @ W @ LambdaInverse
            X = matrixExponentiate(-Lambda * self.source.k0 * self.thickness)
            return (V, W, X)

    def S_matrix(self):
        if self.thickness > 0:
            return self._S_matrix_internal()
        elif self.thickness == 0:

            if self.incident:
                return self._S_matrix_reflection()
            elif self.transmission:
                return self._S_matrix_transmission()
            else:
                raise ValueError('''Semi-infinite film appears to be neither incident or transmissive. 
                Cannot compute S-matrix''')

    def _S_matrix_internal(self):
        (Vi, Wi, Xi) = self.VWX_matrices()
        Ai = A_matrix(Wi, self.Wg, Vi, self.Vg)
        Bi = B_matrix(Wi, self.Wg, Vi, self.Vg)
        Di = D_matrix(Ai, Bi, Xi)

        Si = calculateInternalSMatrixFromRaw(Ai, Bi, Xi, Di);
        return Si;

    def _S_matrix_reflection(self):
        if isinstance(self.Kx, np.ndarray):
            return self._S_matrix_reflection_general()
        else:
            return self._S_matrix_reflection_homogenous()

    def _S_matrix_reflection_homogenous(self):
        (Vi, Wi, X) = self.VWX_matrices()
        Ai = A_matrix(self.Wg, Wi, self.Vg, Vi)
        Bi = B_matrix(self.Wg, Wi, self.Vg, Vi)

        Si = calculateReflectionRegionSMatrixFromRaw(Ai, Bi)
        return Si

    def _S_matrix_reflection_general(self):
        KDimension = self.Kx.shape[0]
        lambdaRef = complexZeros((KDimension*2, KDimension*2))
        Wi = complexIdentity(KDimension * 2)
        Q = self.Q_matrix()

        # I have no idea why we conjugate ur * er and then conjugate the whole thing.
        Kz = conj(sqrt (conj(self.er * self.ur) * \
                        complexIdentity(KDimension) - self.Kx @ self.Kx - self.Ky @ self.Ky))
        lambdaRef[:KDimension, :KDimension] = 1j*Kz
        lambdaRef[KDimension:, KDimension:] = 1j*Kz
        Vi = Q @ pinv(lambdaRef)
        Ai = A_matrix(self.Wg, Wi, self.Vg, Vi)
        Bi = B_matrix(self.Wg, Wi, self.Vg, Vi)

        Sref = calculateReflectionRegionSMatrixFromRaw(Ai, Bi)
        return Sref

    def _S_matrix_transmission(self):
        if isinstance(self.Kx, np.ndarray):
            return self._S_matrix_transmission_general()
        else:
            return self._S_matrix_transmission_homogenous()

    def _S_matrix_transmission_homogenous(self):
        (Vi, Wi, X) = self.VWX_matrices()
        Ai = A_matrix(self.Wg, Wi, self.Vg, Vi);
        Bi = B_matrix(self.Wg, Wi, self.Vg, Vi);

        Si = calculateTransmissionRegionSMatrixFromRaw(Ai, Bi);
        return Si;

    def _S_matrix_transmission_general(self):
        KDimension = self.Kx.shape[0]
        lambdaRef = complexZeros((KDimension*2, KDimension*2))
        Wi = complexIdentity(KDimension * 2)
        Q = self.Q_matrix()

        # I have no idea why we conjugate ur * er and then conjugate the whole thing.
        Kz = conj(sqrt (conj(self.er * self.ur) * complexIdentity(KDimension) - self.Kx @ self.Kx - self.Ky @ self.Ky))
        lambdaRef[:KDimension, :KDimension] = 1j*Kz
        lambdaRef[KDimension:,KDimension:] = 1j*Kz
        Vi = Q @ pinv(lambdaRef)
        Ai = A_matrix(self.Wg, Wi, self.Vg, Vi)
        Bi = B_matrix(self.Wg, Wi, self.Vg, Vi)

        Strn = calculateTransmissionRegionSMatrixFromRaw(Ai, Bi)
        return Strn
