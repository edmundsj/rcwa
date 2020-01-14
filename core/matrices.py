import numpy as np
import scipy as sp
import scipy.linalg

inv = np.linalg.inv;
expm = sp.linalg.expm;
sqrtm = sp.linalg.sqrtm;
sqrt = np.lib.scimath.sqrt; # Takes sqrt of complex numbers
sq = np.square;
eig = sp.linalg.eig # Performs eigendecomposition of identity intuitively (vectors are unit vectors)
norm = np.linalg.norm;
sin = np.sin;
cos = np.cos;

OUTER_BLOCK_SHAPE = (2,2);
PQ_SHAPE = (2,2); # The shape of our core PQ matrices.
TOTAL_SHAPE = OUTER_BLOCK_SHAPE + PQ_SHAPE;
DBGLVL = 2;

# This function makes a couple really subtle assumptions. First, we have to define what TE should be
# when theta = 0. By convention, of how we define theta and phi, this should point along the y-axis.
def aTEM_gen(kx_n, ky_n, kz_n):
    """
    Generates the aTE and aTM vectors from the known kx, ky, kz incident vectors, assuming that
    our device is planar in the x/y direction.
    """
    z_vec_norm = np.array([0,0,-1]); # The z-normal vector (plane wave opposite direction as surface)
    kn_vec = np.array([kx_n, ky_n, kz_n]);
    epsilon = 1e-5; # a small number.

    if(kx_n == 0 and ky_n == 0):
        aTE = np.array([0,1,0]); # This is assuming W is the identity matrix.
    else:
        # HACK. I DO NOT KNOW WHY THE MINUS SIGN IS NEEDED.
        aTE = - np.cross(z_vec_norm, kn_vec);
        aTE = aTE / norm(aTE);

    # I'm pretty sure Rumpf is wrong on this formula. For theta=0, phi=0, aTE should be pointing in
    # the +y direction. If we want aTM to point in the +x direction, and the k-vector is pointing along
    # z, we need to cross aTE with kn, not the other way around.
    aTM = np.cross(aTE, kn_vec);
    aTM = aTM / norm(aTM);

    # Now, we want to shrink down these vectors so that they only contain the x/y polarization
    # information, because we don't use the rest.

    return (aTE[0:2], aTM[0:2]);

def redhefferProduct(SA, SB):
    """ Computes the redheffer star product of
    two matrices A and B. A and B can themselves (I think)
    be matrices. The matrices must be in 2x2 block form for this
    to work properly.
    """

    mat_shape = SA.shape;
    # First, check to make sure SA and SB are the same shape.
    if(mat_shape != SB.shape):
        raise Exception(f'redhefferProduct: SA and SB are not of the same shape. SA is of shape {SA.shape} and SB is of shape {SB.shape}');

    if(mat_shape[0:2] != OUTER_BLOCK_SHAPE):
        raise Exception(f'redhefferProduct: SA is not of block form. SA has shape {SA.shape}');

    block_shape = mat_shape[2:];

    # Making the assumption that the sub-blocks are square.
    ident_block = np.identity(block_shape[0], dtype=np.cdouble);

    SAB = np.zeros(mat_shape, dtype = np.cdouble);
    D = DredGen(SA[0,1], SA[1,1], SB[0,0]);
    F = FredGen(SA[1,1], SB[0,0], SB[1,0]);

    SAB[0,0] = SA[0,0] + D @ SB[0,0] @ SA[1,0];
    SAB[0,1] = D @ SB[0,1];
    SAB[1,0] = F @ SA[1,0];
    SAB[1,1] = SB[1,1] +  F @ SA[1,1] @ SB[0,1];

    return SAB;

def Pi_gen(kx_n, ky_n, eri, uri):
    """
    Computes the P-matrix for the ith layer, given a known relative permeability ui and relative
    permittivity pi. Assumes kx and ky vectors have been normalized to k0. This will have to be changed
    when both kx, xy, ur, and er are tensors.
    """

    P = np.zeros(PQ_SHAPE, dtype=np.cdouble)
    P[0,0] = kx_n*ky_n;
    P[0,1] = uri*eri - np.square(kx_n);
    P[1,0] = sq(ky_n) - uri*eri;
    P[1,1] = - kx_n*ky_n;

    P /= eri;
    return P

def Qi_gen(kx_n, ky_n, eri, uri):
    """
    Computes the Q-matrix for the ith layer, given a known relative permeability ui and relative
    permittivity pi. Assumes kx and ky vectors have been normalized to k0. This will have to be changed
    when both kx, xy, ur, and er are tensors.
    """

    Q = np.zeros(PQ_SHAPE, dtype=np.cdouble)
    Q[0,0] = kx_n * ky_n;
    Q[0,1] = uri*eri - sq(kx_n);
    Q[1,0] = sq(ky_n) - uri*eri;
    Q[1,1] = - kx_n * ky_n;

    Q = Q / uri;
    return Q;

# FUNCTION DOES NOT HAVE UNIT TESTS
def Aij_gen(Wi, Wj, Vi, Vj):
    """
    Computes the matrix Aij for two sets of eigenvector matrices Wi, Wj, Vi, Vj
    """
    return inv(Wi) @ Wj + inv(Vi) @ Vj;

# FUNCTION DOES NOT HAVE UNIT TESTS
def Bij_gen(Wi, Wj, Vi, Vj):
    return inv(Wi) @ Wj - inv(Vi) @ Vj;

def DiGen(Ai, Bi, Xi):
    Ai_inv = inv(Ai);
    return Ai - Xi @ Bi @ Ai_inv @ Xi @ Bi;

def DredGen(S12A, S22A, S11B):
    """
    Generates the D-matrix for the Redheffer star product. NOT the same as the Di matrix.
    """
    return S12A @ inv(np.identity(PQ_SHAPE[0], dtype=np.cdouble) - S11B @ S22A)

def FredGen(S22A, S11B, S21B):
    """
    Generates the F-matrix for computing the Redheffer star product.
    """
    return S21B @ inv(np.identity(PQ_SHAPE[0], dtype=np.cdouble) - S11B @ S22A)

def kzGen(kx_n, ky_n, er, ur):
    return sqrt(er*ur - sq(kx_n) - sq(ky_n));

def Omega_gen(kz_n):
    return np.identity(PQ_SHAPE[0],dtype=np.cdouble) * (0 + 1j)*kz_n;

# FUNCTION DOES NOT HAVE UNIT TESTS. THIS FUNCTION ALSO WILL ONLY WORK FOR LHI MEDIA.
def VWX_gen(kx_n, ky_n, kz_n, er, ur, k0=0, Li=0):
    """
    FUNCTION DOES NOT CURRENTLY WORK. NEEDS TO BE FIXED.
    Generates the V/W matrices (and the X matrix if k0 is nonzero)
    """
    # Unfortunately, for now we have to ignore P and Q, because they are causing numerical
    # instability and severe annoyance. We will instead generate the O2 matrix directly,
    # since it is just a multiple of the identity.
    # Also, eigendecomposition of the identity matrix is an unstable procedure. It will basically
    # give us unit vectors all over the place. For the identity matrix, I will say that we have
    # a single eigenvalue eigval, which is the square root of the multiplier in front of omega squared.
    #P = Pi_gen(kx_n, ky_n, er, ur);
    Q = Qi_gen(kx_n, ky_n, er, ur);

    #O2 = np.ientity(PQ_SHAPE); # Compute the omega squared matrix from the PQ matrix product
    # Unfortunately, eigendecomposition on LH media actually leads to a spread of eigenvectors,
    # because by definition, all possible vectors are vectors of the identity matrix.
    #l2, W = eig(O2); # Eigendecompose this into a diagonal matrix and the basis vectors

    O = Omega_gen(kz_n)

    W = np.identity(PQ_SHAPE[0]); # WRONG WRONG WRONG. THIS IS SUPPOSED TO RELATE X AND Y AND TE/TM MODES.

    O_inv = inv(O);
    X = expm(O * k0 * Li)

    V = Q @ W @ O_inv;

    if(k0 > 0):
        return (V, W, X);
    else:
        return (V, W);

def Si_gen(Ai, Bi, Xi, Di):
    """
    Compute the symmetric scattering matrix using free space (gap layer, Wg)
    The goal is to minimize computation. For each layer, we only want to compute the P/Q/W matrices
    once, and then generate the scattering matrices from that, and return the scattering matrix to
    the main program, which will be used in subsequent computation. The exception is the generation
    of the gap matrices, which we only want to generate once, because they are re-used throughout
    the program
    """
    # The shape of our overall scattering matrix. Will be a matrix of matrices.
    S = np.zeros(TOTAL_SHAPE, dtype=np.cdouble);

    # First, compute all the auxiliary matrices we need to compute our scattering matrix
    Ai_inv = inv(Ai);
    Di_inv = inv(Di);

    S[0,0] = Di_inv @ ((Xi @ Bi @ Ai_inv @ Xi @ Ai) - Bi)
    S[0,1] = Di_inv @ Xi @ (Ai - (Bi @ Ai_inv @ Bi));
    S[1,0] = S[0,1];
    S[1,1] = S[0,0];

    return S;

# NOT COMPLETELY SURE ABOUT THIS FUNCTION. BARELY PASSING UNIT TESTS.
def Sref_gen(Aref, Bref):
    """
    Compute the reflection scattering matrix (the one where we are injecting our excitation wave).
    This matrix is only computed once at the beginning of the simulation.
    """
    S = np.zeros(TOTAL_SHAPE, dtype=np.cdouble);

    Aref_inv = inv(Aref);

    S[0,0] = - Aref_inv @ Bref;
    S[0,1] = 2*Aref_inv;
    S[1,0] = 0.5*(Aref - Bref @ Aref_inv @ Bref)
    S[1,1] = Bref @ Aref_inv;

    return S;

def Strn_gen(Atrn, Btrn):
    """
    Computes the transmission scattering matrix (the one at the 'output' of our device.)
    """

    Atrn_inv = inv(Atrn);
    S = np.zeros(TOTAL_SHAPE, dtype=np.cdouble);

    S[0,0] = Btrn @ Atrn_inv;
    S[0,1] = 0.5* (Atrn - (Btrn @ Atrn_inv @ Btrn))
    S[1,0] = 2* Atrn_inv;
    S[1,1] = - Atrn_inv @ Btrn;

    return S;

def calcEz(kx_n, ky_n, kz_n, Ex, Ey):
    '''
    Calculate the z-component of the electromagnetic field from the x- and y-components using the divergence
    theorem. We are assuming that the material in which we are calculating the z-component is LHI. The Ex
    and Ey components are assumed to be scalars, eri and uri, the relative permittivities and permeabilities,
    are also assumed to be scalars.
    '''
    # First, we calculate kz in our medium from the refractive index and the (ASSUMED) dispersion
    # relation of the medium.
    Ez = (kx_n*Ex + ky_n*Ey) / kz_n
    return Ez;

def calcReflectanceTransmittance(kx_n, ky_n, kz_n, ntrn, ur_ref, ur_trn, Exy_i, Wref, Wtrn, S11, S21):
    '''
    Calculate the reflectance and transmittance given an input electric field vector
    (assumed to be a column array in the form [[Ex],[Ey]]), the incident kz, and the transmitted
    kz.
    '''
    # By default, power conservation is violated.
    R = 0;
    T = 0;

    # First, calculate kz_n given an assumed LHI dispersion relation for the incident and outgoing media.
    kz_ntrn = sqrt(sq(ntrn) - sq(kx_n) - sq(ky_n));

    # Next, calculate the transmitted and reflected fields given the incident fields
    Exy_ref = Wref @ S11 @ inv(Wref) @ Exy_i;
    Exy_trn = Wtrn @ S21 @ inv(Wtrn) @ Exy_i;

    # Now that we have the x and y fields for each region, we can calculate the z-component.
    Ez_i = calcEz(kx_n, ky_n, kz_n, Exy_i[0,0], Exy_i[1,0]);
    Ez_ref = calcEz(kx_n, ky_n, kz_n, Exy_ref[0,0], Exy_ref[1,0]);
    Ez_trn = calcEz(kx_n, ky_n, kz_ntrn, Exy_trn[0,0], Exy_trn[1,0]);

    # Now compose everything into a total 3-dimensional E-field vector.
    Etot_i = np.append(Exy_i, Ez_i);
    Etot_ref = np.append(Exy_ref, Ez_ref);
    Etot_trn = np.append(Exy_trn, Ez_trn);

    # Finally, we can compute the ratio of the E-fields and the real part of kz.
    R = sq(norm(Etot_ref)) / sq(norm(Etot_i));
    T = sq(norm(Etot_trn)) / sq(norm(Etot_i)) * np.real(ur_ref * kz_ntrn / (ur_trn * kz_n));

    if(R + T < 1 - 1e-10):
        print(f"R + T = {R + T}. Do you have an absorbing medium? Power conservation appears to be violated.");
    elif(R + T > 1 + 1e-10):
        print(f"R + T = {R + T}. Do you have a gain medium? Power conservation appears to be violated.");

    return (R, T);

