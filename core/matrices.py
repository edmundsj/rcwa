import numpy as np
import scipy as sp
import scipy.linalg

inv = np.linalg.inv;
matrixExponentiate = sp.linalg.expm
matrixSquareRoot = sp.linalg.sqrtm
sqrt = np.lib.scimath.sqrt; # Takes sqrt of complex numbers
sq = np.square;
eig = sp.linalg.eig # Performs eigendecomposition of identity intuitively (vectors are unit vectors)
norm = np.linalg.norm;
sin = np.sin;
cos = np.cos;

OUTER_BLOCK_SHAPE = (2,2);
scatteringElementShape = (2,2); # The shape of our core PQ matrices.
scatteringMatrixShape = OUTER_BLOCK_SHAPE + scatteringElementShape;
scatteringElementShape = (2,2);
scatteringElementSize = scatteringElementShape[0];
DBGLVL = 2;

class _const:
    class ConstError(TypeError): pass
    class ConstCaseError(ConstError): pass

    def __setattr__(self, name, value):
        if name in self.__dict__:
            raise self.ConstError(f"Can't change const {name}")
        self.__dict__[name] = value;

# This function makes a couple really subtle assumptions. First, we have to define what TE should be
# when theta = 0. By convention, of how we define theta and phi, this should point along the y-axis.
def aTEMGen(kx, ky, kz):
    """
    Generates the aTE and aTM vectors from the known kx, ky, kz incident vectors, assuming that
    our device is planar in the x/y direction.
    """
    deviceNormalUnitVector = complexArray([0,0,-1]); # The z-normal vector (plane wave opposite direction as surface)
    kn_vec = complexArray([kx, ky, kz]);
    epsilon = 1e-3; # a small number. If both our kx and ky vectors are tiny, our cross product will not
    # be computed properly, and we need to fix that.

    if(abs(kx) < epsilon and abs(ky) < epsilon):
        aTE = np.array([0,1,0]); # This is assuming W is the identity matrix.
    else:
        # HACK. I DO NOT KNOW WHY THE MINUS SIGN IS NEEDED.
        aTE = - np.cross(deviceNormalUnitVector, kn_vec);
        aTE = aTE / norm(aTE);

    # I'm pretty sure Rumpf is wrong on this formula. For theta=0, phi=0, aTE should be pointing in
    # the +y direction. If we want aTM to point in the +x direction, and the k-vector is pointing along
    # z, we need to cross aTE with kn, not the other way around.
    aTM = np.cross(aTE, kn_vec);
    aTM = aTM / norm(aTM);

    # Now, we want to shrink down these vectors so that they only contain the x/y polarization
    # information, because we don't use the rest.

    return (aTE[0:2], aTM[0:2]);

def complexArray(arrayInListForm):
    """ Wrapper for numpy array declaration that forces arrays to be complex doubles """
    return np.array(arrayInListForm, dtype=np.cdouble);

def complexIdentity(matrixSize):
    """ Wrapper for numpy identity declaration that forces arrays to be complex doubles """
    return np.identity(matrixSize, dtype=np.cdouble);

def complexZeros(matrixDimensionsTuple):
    """ Wrapper for numpy zeros declaration that forces arrays to be complex doubles """
    return np.zeros(matrixDimensionsTuple, dtype=np.cdouble);

def generateTransparentSMatrix():
    STransparent = complexZeros(scatteringMatrixShape);
    STransparent[0,1] = complexIdentity(scatteringElementSize);
    STransparent[1,0] = complexIdentity(scatteringElementSize);
    return STransparent;

def calculateRedhefferProduct(SA, SB):
    """
    Computes the redheffer star product of
    two matrices A and B. A and B can themselves (I think)
    be matrices. The matrices must be in 2x2 block form for this
    to work properly.
    """

    mat_shape = SA.shape;
    # First, check to make sure SA and SB are the same shape.
    if(mat_shape != SB.shape):
        raise Exception(f'redhefferProduct: SA and SB are not of the same shape. SA is of shape {SA.shape} and SB is of shape {SB.shape}');

    # Making the assumption that the sub-blocks are square.
    SAB = complexZeros(mat_shape);
    D = calculateRedhefferDMatrix(SA[0,1], SA[1,1], SB[0,0]);
    F = calculateRedhefferFMatrix(SA[1,1], SB[0,0], SB[1,0]);

    SAB[0,0] = SA[0,0] + D @ SB[0,0] @ SA[1,0];
    SAB[0,1] = D @ SB[0,1];
    SAB[1,0] = F @ SA[1,0];
    SAB[1,1] = SB[1,1] +  F @ SA[1,1] @ SB[0,1];

    return SAB;

def calculatePMatrix(kx, ky, eri, uri):
    """
    Computes the P-matrix for the ith layer, given a known relative permeability ui and relative
    permittivity pi. Assumes kx and ky vectors have been normalized to k0. This will have to be changed
    when both kx, xy, ur, and er are tensors.
    """
    P = complexZeros(scatteringElementShape);
    P[0,0] = kx*ky;
    P[0,1] = uri*eri - np.square(kx);
    P[1,0] = sq(ky) - uri*eri;
    P[1,1] = - kx*ky;

    P /= eri;
    return P

def calculateQMatrix(kx, ky, eri, uri):
    """
    Computes the Q-matrix for the ith layer, given a known relative permeability ui and relative
    permittivity pi. Assumes kx and ky vectors have been normalized to k0. This will have to be changed
    when both kx, xy, ur, and er are tensors.
    """

    Q = complexZeros(scatteringElementShape);
    Q[0,0] = kx * ky;
    Q[0,1] = uri*eri - sq(kx);
    Q[1,0] = sq(ky) - uri*eri;
    Q[1,1] = - kx * ky;

    Q = Q / uri;
    return Q;

# FUNCTION DOES NOT HAVE UNIT TESTS
def calculateScatteringAMatrix(Wi, Wj, Vi, Vj):
    """
    Computes the matrix Aij for two sets of eigenvector matrices Wi, Wj, Vi, Vj
    """
    return inv(Wi) @ Wj + inv(Vi) @ Vj;

def calculateScatteringBMatrix(Wi, Wj, Vi, Vj): # UNIT TESTS COMPLETE
    return inv(Wi) @ Wj - inv(Vi) @ Vj;

def calculateScatteringDMatrix(Ai, Bi, Xi): # UNIT TESTS COMPLETE
    AiInverse = inv(Ai);
    return Ai - Xi @ Bi @ AiInverse @ Xi @ Bi;

def calculateRedhefferDMatrix(S12A, S22A, S11B): # UNIT TESTS COMPLETE
    """
    Generates the D-matrix for the Redheffer star product. NOT the same as the Di matrix.
    """
    return S12A @ inv(complexIdentity(scatteringElementShape[0]) - S11B @ S22A)

def calculateRedhefferFMatrix(S22A, S11B, S21B): # UNIT TESTS COMPLETE
    """
    Generates the F-matrix for computing the Redheffer star product.
    """
    return S21B @ inv(complexIdentity(scatteringElementShape[0]) - S11B @ S22A)

def calculateKz(kx, ky, er, ur): # UNIT TESTS COMPLETE
    return sqrt(er*ur - sq(kx) - sq(ky));

def calculateKVector(theta, phi, er, ur):
    n = sqrt(er*ur);
    kx = n * sin(theta) * cos(phi);
    ky = n * sin(theta) * sin(phi);
    kz = n * cos(theta);
    return complexArray([kx, ky, kz]);

def calculateOmegaMatrix(kz): # UNIT TESTS COMPLETE
    return complexIdentity(2)* (0 + 1j)*kz;

def calculateVWXMatrices(kx, ky, kz, er, ur, k0=0, Li=0): # UNIT TESTS COMPLETE
    """
    FUNCTION DOES NOT CURRENTLY WORK. NEEDS TO BE FIXED.
    Generates the V/W matrices (and the X matrix if k0 is nonzero)
    """
    Q = calculateQMatrix(kx, ky, er, ur);

    O = calculateOmegaMatrix(kz);

    W = complexIdentity(scatteringElementShape[0]);

    OInverse = inv(O);
    X = matrixExponentiate(O * k0 * Li)

    V = Q @ W @ OInverse;

    if(k0 > 0):
        return (V, W, X);
    else:
        return (V, W);


def calculateInternalSMatrix(kx, ky, er, ur, k0, Li, Wg, Vg):

    # First, calculate the kz component inside this layer
    kz = calculateKz(kx, ky, er, ur);
    (Vi, Wi, Xi) = calculateVWXMatrices(kx, ky, kz, er, ur, k0, Li);

    Ai = calculateScatteringAMatrix(Wi, Wg, Vi, Vg);
    Bi = calculateScatteringBMatrix(Wi, Wg, Vi, Vg);

    Di = calculateScatteringDMatrix(Ai, Bi, Xi);

    Si = calculateInternalSMatrixFromRaw(Ai, Bi, Xi, Di);
    return Si;

def calculateReflectionRegionSMatrix(kx, ky, er, ur, Wg, Vg):
    """
    Calculates S-matrix for reflection region using raw material data
    """
    kz = calculateKz(kx, ky, er, ur);
    (Vi, Wi) = calculateVWXMatrices(kx, ky, kz, er, ur);

    Ai = calculateScatteringAMatrix(Wg, Wi, Vg, Vi);
    Bi = calculateScatteringBMatrix(Wg, Wi, Vg, Vi);

    Si = calculateReflectionRegionSMatrixFromRaw(Ai, Bi);
    return Si;

def calculateTransmissionRegionSMatrix(kx, ky, er, ur, Wg, Vg):
    """
    Calculates S-matrix for reflection region using raw material data
    """
    kz = calculateKz(kx, ky, er, ur);
    (Vi, Wi) = calculateVWXMatrices(kx, ky, kz, er, ur);

    Ai = calculateScatteringAMatrix(Wg, Wi, Vg, Vi);
    Bi = calculateScatteringBMatrix(Wg, Wi, Vg, Vi);

    Si = calculateTransmissionRegionSMatrixFromRaw(Ai, Bi);
    return Si;

def calculateInternalSMatrixFromRaw(Ai, Bi, Xi, Di):
    """
    Compute the symmetric scattering matrix using free space (gap layer, Wg)
    The goal is to minimize computation. For each layer, we only want to compute the P/Q/W matrices
    once, and then generate the scattering matrices from that, and return the scattering matrix to
    the main program, which will be used in subsequent computation. The exception is the generation
    of the gap matrices, which we only want to generate once, because they are re-used throughout
    the program
    """
    # The shape of our overall scattering matrix. Will be a matrix of matrices.
    S = complexZeros(scatteringMatrixShape);

    # First, compute all the auxiliary matrices we need to compute our scattering matrix
    AiInverse = inv(Ai);
    DiInverse = inv(Di);

    S[0,0] = DiInverse @ ((Xi @ Bi @ AiInverse @ Xi @ Ai) - Bi)
    S[0,1] = DiInverse @ Xi @ (Ai - (Bi @ AiInverse @ Bi));
    S[1,0] = S[0,1];
    S[1,1] = S[0,0];

    return S;

def calculateReflectionRegionSMatrixFromRaw(AReflectionRegion, BReflectionRegion):
    """
    """
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
    """
    Computes the transmission scattering matrix (the one at the 'output' of our device.) from the raw
    """
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
    '''
    Calculate the z-component of the electromagnetic field from the x- and y-components using the divergence
    theorem. We are assuming that the material in which we are calculating the z-component is LHI. The Ex
    and Ey components are assumed to be scalars, eri and uri, the relative permittivities and permeabilities,
    are also assumed to be scalars.
    '''
    # First, we calculate kz in our medium from the refractive index and the (ASSUMED) dispersion
    # relation of the medium.
    Ez = - (kx*Ex + ky*Ey) / kz
    return Ez;

def calculateRT(kzReflectionRegion, kzTransmissionRegion,
        urReflectionRegion, urTransmissionRegion, ExyzReflected, ExyzTransmitted):
    '''
    Calculate the reflectance and transmittance given an input electric field vector
    (assumed to be a column array in the form [[Ex],[Ey]]), the incident kz, and the transmitted
    kz. WARNING: assumes incident fields are normalized to a magnitude of 1. We need to enforce
    this elsewhere.
    '''
    R = sq(norm(ExyzReflected))
    T = sq(norm(ExyzTransmitted))*np.real(kzTransmissionRegion / urTransmissionRegion) / \
            (kzReflectionRegion / urReflectionRegion);

    return (R, T);

