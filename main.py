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

# First, figure out how many layers we have
num_internal_layers = len(t);
if(len(er) != num_internal_layers + 2):
    raise Exception(f"Error: The number of layers is not equal to the number of permittivities. Number of permittivities is {len(er)} and number of layers is {num_layers+2}");

if(len(ur) != num_internal_layers + 2):
    raise Exception(f"Error: The number of layers is not equal to the number of permeabilities. Number of permeabilities is {len(ur)} and number of layers is {num_layers+2}");


print("Initializing Simulation... Setting up materials, polarization, incident wave...")
# Setup material parameters used in the simulation
urReflectionRegion = ur[0];
urTransmissionRegion = ur[-1];
erReflectionRegion = er[0];
erTransmissionRegion = er[-1];
nReflectionRegion = sqrt(erReflectionRegion * urReflectionRegion);
nTransmissionRegion = sqrt(erReflectionRegion * urReflectionRegion);

# Setup wavelength, k vector, polarization state, and polarization vector
wavelength = sources[0][0]; # The wavelength of our simulation
k0 = 2*np.pi / wavelength;
theta = sources[0][1]; # For now, we will only have one source.
phi = sources[0][2];

kx = nReflectionRegion * sin(theta)*cos(phi);
ky = nReflectionRegion * sin(theta)*sin(phi);
kzReflectionRegion = nReflectionRegion * cos(theta);
kzTransmissionRegion = sqrt(erTransmissionRegion * urTransmissionRegion - sq(kx) - sq(ky));
aTE, aTM = aTEMGen(kx, ky, kzReflectionRegion);

pTEM = complexArray([sources[0][3], sources[0][4]]);
pTEM = pTEM / norm(pTEM); # enforce normalized polarization state
(pTE, pTM) = (pTEM[0], pTEM[1]);
ExyIncident = pTE*aTE + pTM*aTM;

# Generate our gap matrices once so we don't have to keep re-generating them.
erg = 1 + sq(kx) + sq(ky); # Our gap permittivity (just the permittivity of free space)
urg = 1; # Our gap permeability (just the permeability of free space)
VGap, WGap = calculateVWXMatrices(kx, ky, kzReflectionRegion, erg, urg);

# Generate our reference and transmission scattering matrices 
SReflectionRegion = calculateReflectionRegionSMatrix(kx, ky,
        erReflectionRegion, urReflectionRegion, WGap, VGap);

STransmissionRegion = calculateTransmissionRegionSMatrix(kx, ky,
        erTransmissionRegion, urTransmissionRegion, WGap, VGap);

# Initialize our system's scattering matrix. We should have no reflection in either way and total transmission
SGlobal = generateTransparentSMatrix();

if(DBGLVL >= 2):
    print(f""""Initial Variables --- \nSGlobal: {SGlobal}\n\nSref:\n{SReflectionRegion}\n\n
            STransmissionRegion:\n{STransmissionRegion}\n""")

print("Generating scattering matrices for internal layers, updating system matrix...");
for i in range(num_internal_layers):
    # Generate the ith scattering matrix for layer i
    Si = calculateInternalSMatrix(kx, ky, er[i+1], ur[i+1], k0, t[i], WGap, VGap);

    if(DBGLVL >= 2):
        print(f"LAYER {i} ---\nS{i}:\n{Si}\n\n");

    # Update our system's scattering matrix using the redheffer star product with the ith scattering matrix
    SGlobal= calculateRedhefferProduct(SGlobal, Si);
    print(f"Layer {i+1} Done");

# Update the global scattering matrix to include the reflection and transmission regions
SGlobal = calculateRedhefferProduct(SReflectionRegion, SGlobal);
SGlobal = calculateRedhefferProduct(SGlobal, STransmissionRegion);

# Now, we need to calculate the reflected and transmitted fields using our S-matrices
ExyReflected = SGlobal[0,0] @ ExyIncident;
ExyTransmitted = SGlobal[1,0] @ ExyIncident;
EzReflected = calculateEz(kx, ky, kzReflectionRegion, ExyReflected[0], ExyReflected[1]);
EzTransmitted = calculateEz(kx, ky, kzTransmissionRegion, ExyTransmitted[0], ExyTransmitted[1]);
ExyzReflected = np.append(ExyReflected, EzReflected);
ExyzTransmitted = np.append(ExyTransmitted, EzTransmitted);

if(DBGLVL >= 2):
    print(f"ExyzIncident: {ExyIncident}\nExyReflected: {ExyzReflected}\nExyzTransmitted: {ExyzTransmitted}");

# And at long last, calculate the reflectance and transmittance of our device.
R, T = calculateRT(kzReflectionRegion, kzTransmissionRegion, urReflectionRegion, urTransmissionRegion,
        ExyzReflected, ExyzTransmitted);

print("Done with all layers");
print(f"R: {R}, T: {T}, R+T={R+T}");
if(DBGLVL >= 2):
    print(f"\nSGlobal:\n{SGlobal}\n");
