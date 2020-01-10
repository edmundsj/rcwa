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
DBGLVL = 2;

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
    D = SA[0,1] @ inv(ident_block - SB[0,0] @ SA[1,1]);
    F = SB[1,0] @ inv(ident_block - SA[1,1] @ SB[0,0]);

    SAB[0,0] = SA[0,0] + D @ SB[0,0] @ SA[1,0];
    SAB[0,1] = D @ SB[0,1];
    SAB[1,0] = F @ SA[1,0];
    SAB[1,1] = SB[1,1] + SB[1,0] @ F @ SA[1,1] @ SB[0,1];

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
    Q[0,1] = uri*eri - np.square(kx_n);
    Q[1,0] = np.square(ky_n) - uri*eri;
    Q[1,1] = - kx_n * ky_n;

    Q /= uri;
    return Q;

# Computes the matrix Aij for two sets of eigenvector matrices Wi, Wj, Vi, Vj
def Aij_gen(Wi, Wj, Vi, Vj):
    return inv(Wi) @ Wj + inv(Vi) @ Vj;

def Bij_gen(Wi, Wj, Vi, Vj):
    return inv(Wi) @ Wj - inv(Vi) @ Vj;

def Xi_gen(lambda_i, k0, Li):
    return expm(lambda_i*k0*Li);

def Si_gen(k0, Li, kx_n, ky_n, eri, uri, Wg, Vg):
    """
    Compute the symmetric scattering matrix using free space (gap layer, Wg)
    The goal is to minimize computation. For each layer, we only want to compute the P/Q/W matrices
    once, and then generate the scattering matrices from that, and return the scattering matrix to
    the main program, which will be used in subsequent computation. The exception is the generation
    of the gap matrices, which we only want to generate once, because they are re-used throughout
    the program
    """
    Pi = Pi_gen(kx_n, ky_n, eri, uri);
    Qi = Qi_gen(kx_n, ky_n, eri, uri);

    O2i = Pi @ Qi; # Omega-squared matrix that we want to take the matrix exponential of
    l2i, Wi = eig(O2i); # Get W- and lambda-squared matrices from the eigendecomposition
    lambda_i = np.diag(sqrt(l2i));
    Vi = Qi @ Wi @ inv(lambda_i); #Compute magnetic field eigenvectors from definition

    inner_shape = Pi.shape;
    total_shape = OUTER_BLOCK_SHAPE + inner_shape;

    # The shape of our overall scattering matrix. Will be a matrix of matrices.
    S = np.zeros(total_shape, dtype=np.cdouble);

    # First, compute all the auxiliary matrices we need to compute our scattering matrix
    Ai = Aij_gen(Wi, Wg, Vi, Vg);
    Bi = Bij_gen(Wi, Wg, Vi, Vg);
    Xi = Xi_gen(lambda_i, k0, Li);
    Ai_inv = inv(Ai);
    block_1 = inv(Ai - Xi @ Bi @ Ai_inv @ Xi @ Bi);

    if(DBGLVL >= 2):
        print("Calling Si_gen():")
        print(f"Pi:\n{Pi}\nQi:\n{Qi}\nlambda_i:\n{lambda_i}\nVi:\n{Vi}\nAi:{Ai}\n\nBi:\n{Bi}\n\nXi:\n{Xi}")

    S[0,0] = block_1 @ (Xi @ Bi @ Ai_inv @ Xi @ Ai - Bi)
    S[0,1] = block_1 @ Xi @ (Ai - Bi @ Ai_inv @ Bi);
    S[1,0] = S[0,1];
    S[1,1] = S[0,0];

    return S;

def SWref_gen(kx_n, ky_n, er_ref, ur_ref, Wg, Vg):
    """
    Compute the reflection scattering matrix (the one where we are injecting our excitation wave).
    This matrix is only computed once at the beginning of the simulation.
    """
    inner_shape = PQ_SHAPE;
    total_shape = OUTER_BLOCK_SHAPE + inner_shape;

    Pref = Pi_gen(kx_n, ky_n, er_ref, ur_ref);
    Qref = Qi_gen(kx_n, ky_n, er_ref, ur_ref);
    O2ref = Pref @ Qref;
    l2ref, Wref = eig(O2ref);

    lambda_ref = np.diag(sqrt(l2ref));

    Vref = Qref @ Wref @ inv(lambda_ref);

    S = np.zeros(total_shape, dtype=np.cdouble);

    # I don't trust this line of code (from slide 7 secture 2c of CEM lecture notes)
    # This is inconsistent with prior notation of Aij and Bij, but I'm going to write it down as-is
    Aref = Aij_gen(Wg, Wref, Vg, Vref);
    Bref = Bij_gen(Wg, Wref, Vg, Vref);

    Aref_inv = inv(Aref);

    S[0,0] = - Aref_inv @ Bref;
    S[0,1] = 2*Aref_inv;
    S[1,0] = 0.5*(Aref - Bref @ Aref_inv @ Bref)
    S[1,1] = Bref @ Aref_inv;

    return (S, Wref);

def SWtrn_gen(kx_n, ky_n, er_trn, ur_trn, Wg, Vg):
    """
    Computes the transmission scattering matrix (the one at the 'output' of our device.)
    """
    Ptrn = Pi_gen(kx_n, ky_n, er_trn, ur_trn);
    Qtrn = Qi_gen(kx_n, ky_n, er_trn, ur_trn);
    O2trn = Ptrn @ Qtrn;
    l2trn, Wtrn = eig(O2trn);

    lambda_trn = np.diag(sqrt(l2trn));

    Vtrn = Qtrn @ Wtrn @ inv(lambda_trn);

    inner_shape = Wtrn.shape;
    total_shape = OUTER_BLOCK_SHAPE + inner_shape;

    Atrn = Aij_gen(Wg, Wtrn, Vg, Vtrn);
    Btrn = Bij_gen(Wg, Wtrn, Vg, Vtrn);
    Atrn_inv = inv(Atrn);

    S = np.zeros(total_shape, dtype=np.cdouble);
    S[0,0] = Btrn @ Atrn_inv;
    S[0,1] = 0.5* (Atrn - Btrn @ Atrn_inv @ Btrn)
    S[1,0] = 2* Atrn_inv;
    S[1,1] = - Atrn_inv @ Btrn;

    return (S, Wtrn);

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

