# Author: Jordan Edmunds, Ph.D. Student, UC Berkeley
# Contact: jordan.e@berkeley.edu
# Creation Date: 11/01/2019
#
# TODO:
# Fix netlist parser so it can handle zero layers and assume the input and output layers are free space
# (Eventually, maybe never):
# - relax the constraint that the transmitting and incident medium must be LHI materials.

import numpy as np
import scipy as sp
import scipy.linalg
import sys
from core.matrices import *
from netlist.netlist_parser import *
import matplotlib.pyplot as plt

# 1. The class NetlistParser parses a netlist and turns everything into a "Mask" or "Field" object.
# The masks and field are returned so that they are sorted in ascending order with
# respect to their coordinate on the optical axis.
arguments = len(sys.argv) - 1; # The number of arguments
netlist_location = './netlist/sample_netlist.txt';
print(arguments);
print(sys.argv);
if(arguments > 0):
    print(f"Using user defined netlist {sys.argv[1]}")
    netlist_location = sys.argv[1];
print("Parsing netlist... ");
parser = NetlistParser(netlist_location);

# Up until this point, everything is totally general, but now we have to decide
# what it is our parser returns. I will try to have a consistent netlist format across
# all my codes to the degree that is possible.
[er, ur, t, sources] = parser.parseNetlist();
print(f"Done. Found:\nPermittivities:{er}\nPermeabilities:{ur}\nInternal Layers: {t}")

# Now that we have extracted our permittivities, permeabilities, thicknesses, and sources,
# we are ready to start a simulation. We do that first by initializing our scattering matrix.
# For now, we formulated everything using 2x2 matrices, so we will put things in block matrix form.
# NOTE: HOW NUMPY DOES MULTIPLICATION MEANS WE CANNOT DIRECTLY MULTIPLY THESE HIGHER_DIMENSIONAL MATRICES.

# First, figure out how many layers we have
num_internal_layers = len(t);
if(len(er) != num_internal_layers + 2):
    raise Exception(f"Error: The number of layers is not equal to the number of permittivities. Number of permittivities is {len(er)} and number of layers is {num_layers+2}");

if(len(ur) != num_internal_layers + 2):
    raise Exception(f"Error: The number of layers is not equal to the number of permeabilities. Number of permeabilities is {len(ur)} and number of layers is {num_layers+2}");


# Now, figure out from phi and theta what our incident wavevector is, normalized to the wavevector of free
# space
ur_ref = ur[0];
ur_trn = ur[-1];
er_ref = er[0];
er_trn = er[-1];
nref = sqrt(er_ref * ur_ref); # Get the refractive index for our first semi-infinite layer.
ntrn = sqrt(er_trn * ur_trn); # Refractive index for our final semi-infinite layer

wavelength = sources[0][0]; # The wavelength of our simulation
k0 = 2*np.pi / wavelength;
theta = sources[0][1]; # For now, we will only have one source.
phi = sources[0][2];

# The normalized directional cosines multipled by the refractive index of the incident material
kx_n = nref * np.sin(theta)*np.cos(phi);
ky_n = nref * np.sin(theta)*np.sin(phi);
kz_n = nref * np.cos(theta);
kn = np.array([kx_n, ky_n, kz_n]);
z_vec = np.array([0,0,1]);

# Calculate the TE / TM vectors, our polarization state, and the xy electric field.
if(theta == 0):
    aTE = np.array([0,1,0]);
else:
    aTE = np.cross(kn, z_vec) / norm(np.cross(kn, z_vec)); # the TE field vector (normalized)

aTM = np.cross(aTE, kn) / norm(np.cross(aTE, kn)); # the TM field vector (normalized)
pTE = sources[0][3];
pTM = sources[0][4];
Exy_i = pTE*aTE + pTM*aTM;
Exy_i = Exy_i[0:2];
Exy_i = Exy_i.reshape(2,1);

# Generate our gap matrices once so we don't have to keep re-generating them.
print("Initializing Simulation... Generating gap, transmission, reflection, system scattering matrices...");
erg = 1; # Our gap permittivity (just the permittivity of free space)
urg = 1; # Our gap permeability (just the permeability of free space)

Pg = Pi_gen(kx_n, ky_n, erg, urg);
Qg = Qi_gen(kx_n, ky_n, erg, urg);
O2g = Pg @ Qg;
l2g, Wg = eig(O2g);

lambda_g = np.diag(sqrt(l2g));
Vg = Qg @ Wg @ inv(lambda_g);

# Generate our reference and transmission scattering matrices 
(Sref, Wref) = SWref_gen(kx_n, ky_n, er_ref, ur_ref, Wg, Vg);
(Strn, Wtrn) = SWtrn_gen(kx_n, ky_n, er_trn, ur_trn, Wg, Vg);

# Initialize our system's scattering matrix. We should have no reflection in either way and total transmission
total_shape = OUTER_BLOCK_SHAPE + PQ_SHAPE
Sglobal = np.zeros(total_shape, dtype = np.cdouble);
Sglobal[1, 0] = np.identity(PQ_SHAPE[0], dtype=np.cdouble); # Set our S21 coefficient to 1 (total transmission)
Sglobal[0,1] = np.identity(PQ_SHAPE[0], dtype=np.cdouble); # Set our S12 coefficient to 1 (total backwards transmission)

Sdevice = np.zeros(total_shape, dtype = np.cdouble);
Sdevice[1, 0] = np.identity(PQ_SHAPE[0], dtype=np.cdouble); # Set our S21 coefficient to 1 (total transmission)
Sdevice[0,1] = np.identity(PQ_SHAPE[0], dtype=np.cdouble); # Set our S12 coefficient to 1 (total backwards transmission)

if(DBGLVL >= 2):
    print(f"Initial Variables --- \nSglobal: {Sglobal}\n\nSdevice: \n{Sdevice}\nSref:\n{Sref}\n\nStrn:\n{Strn}\n")

print("Generating scattering matrices for internal layers, updating system matrix...");
for i in range(num_internal_layers):
    # Generate the ith scattering matrix for layer i
    Si = Si_gen(k0, t[i], kx_n, ky_n, er[i+1], ur[i+1], Wg, Vg);

    if(DBGLVL >= 2):
        print(f"LAYER {i} ---\nS{i}:\n{Si}\n\n");

    # Update our system's scattering matrix using the redheffer star product with the ith scattering matrix
    Sdevice = redhefferProduct(Sdevice, Si);
    print(f"Layer {i+1} Done");

# Finally, find the total scattering matrix.
Sglobal = redhefferProduct(Sref, Sdevice);
Sglobal = redhefferProduct(Sglobal, Strn);

# And at long last, calculate the reflectance and transmittance of our device.
[R, T] = calcReflectanceTransmittance(kx_n, ky_n, kz_n, ntrn, ur_ref, ur_trn, Exy_i, Wref, Wtrn, \
        Sglobal[0,0], Sglobal[1,0]);
print("Done with all layers");
print(f"R: {R}, T: {T}, R+T={R+T}");
if(DBGLVL >= 2):
    print(f"Sdevice:\n{Sdevice}\n\nSglobal:\n{Sglobal}\n");
