# Fresnel equations definitions for TE and TM modes
import context
from matrices import *

def rTE(kz1, kz2, er1, er2, ur1, ur2):
    return (ur2 * kz1 - ur1 * kz2) / (ur2 * kz1 + ur1 * kz2);

# This appears to be correct
def rTM(kz1, kz2, er1, er2, ur1, ur2):
    return (er1 * kz2 - er2 * kz1) / (er1 * kz2 + er2 * kz1);

def tTE(kz1, kz2, er1, er2, ur1, ur2):
    return 1 + rTE(kz1, kz2, er1, er2, ur1, ur2);

def tTM(kz1, kz2, er1, er2, ur1, ur2):
    eta_1 = sqrt(ur1 / er1);
    eta_2 = sqrt(ur2 / er2);
    n1 = sqrt(ur1 * er1);
    n2 = sqrt(ur2 * er2);
    return eta_2 / eta_1 * (1 - rTM(kz1, kz2, er1, er2, ur1, ur2));
    #return kz1 / kz2 * n2 / n1 * (1 + rTM(kz1, kz2, er1, er2, ur1, ur2));
    # return n1 / n2 * (1 + rTM(kz1, kz2, er1, er2, ur1, ur2));

def fresnelSMatrixInterface(kx_n, ky_n, er1, er2, ur1, ur2):
    n1_sq = er1 * ur1;
    n2_sq = er2 * ur2;

    # Normalized ki parameters. k1t stands for transverse.
    kt_sq = sq(kx_n) + sq(ky_n); # Transverse k-vector. Will be the same in all media.
    kz1 = sqrt(n1_sq - kt_sq); # z-component of k vector in media 1
    kz2 = sqrt(n2_sq - kt_sq); # z-component of k vector in media 2

    # Amplitude reflection coefficients going from material 1 to material 2
    # This may be a point of contention. This set of equations takes the H-field to be
    # flipped at the interface (so our boundary EQS for TM modes are Hi - Hr = Ht
    # and Ei*cos(th) + Er*cos(th) = Et*cos(th2). We need to make sure this is
    # the same convention as our scattering matrix, otherwise things will be wrong.
    rTE12 = rTE(kz1, kz2, er1, er2, ur1, ur2);
    tTE12 = tTE(kz1, kz2, er1, er2, ur1, ur2);
    rTM12 = rTM(kz1, kz2, er1, er2, ur1, ur2);
    tTM12 = tTM(kz1, kz2, er1, er2, ur1, ur2);

    # Now, we want to compute the amplitude coefficients going back the other way.
    # I believe all we have to do is change our definition of 1 and 2, so we just flip
    # all of the arguments.
    rTE21 = rTE(kz2, kz1, er2, er1, ur2, ur1);
    tTE21 = tTE(kz2, kz1, er2, er1, ur2, ur1);
    rTM21 = rTM(kz2, kz1, er2, er1, ur2, ur1);
    tTM21 = tTM(kz2, kz1, er2, er1, ur2, ur1);

    # Now we are ready to compute our s-parameters matrices in the TE/TM basis.
    # There is no cross-coupling between the TE and TM modes, so the off-diagonal
    # terms of all these matrices are zero.

    # S11 matrix: reflection from material 1 -> material 1
    S11 = np.zeros((2,2), dtype=np.cdouble)
    S11[0,0] = rTM12; # the TE reflection coefficient
    S11[1,1] = rTE12; # the TM reflection coefficient

    # S21 matrix: transmission from material 1 -> material 2
    S21 = np.zeros((2,2), dtype=np.cdouble)
    S21[0,0] = tTM12;
    S21[1,1] = tTE12;

    # S22 matrix: transmission from material 2 -> material 2
    S22 = np.zeros((2,2), dtype=np.cdouble)
    S22[0,0] = rTM21;
    S22[1,1] = rTE21;

    # S12 matrix: transmission from material 2 -> material 1
    S12 = np.zeros((2,2), dtype=np.cdouble)
    S12[0,0] = tTM21;
    S12[1,1] = tTE21;

    return (S11, S12, S21, S22);
