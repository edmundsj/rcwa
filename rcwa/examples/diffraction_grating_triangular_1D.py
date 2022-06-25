import unittest
from rcwa import Source, Layer, LayerStack, Crystal, Solver, TriangularGrating
from rcwa.utils import Plotter
from rcwa.shorthand import complexArray
import numpy as np

def solve_system():
    permittivity_data = np.array([1, 1, 1, 1, 3, 3, 3, 3])
    permeability_data = 1 + 0 * permittivity_data

    reflection_layer = Layer(er=1.0, ur=1.0)
    transmission_layer = Layer(er=9.0, ur=1.0)

    wavelength = 0.5
    deg = np.pi / 180
    k0 = 2*np.pi/wavelength
    theta = 60 * deg
    phi = 1*deg
    pTEM = 1/np.sqrt(2)*complexArray([1,1j])
    source = Source(wavelength=wavelength, theta=theta, phi=phi, pTEM=pTEM, layer=reflection_layer)
    lattice_vector = [1.0, 0, 0]

    crystal_thickness = 0.5

    N_harmonics = 11

    grating = TriangularGrating(period=2, t=0.5, n=4, n_void=1,Nx=500)
    layer_stack = LayerStack(reflection_layer, *grating.slice(), transmission_layer)

    solver_1d = Solver(layer_stack, source, N_harmonics)
    solver_1d.solve()
    return solver_1d

if __name__ == 'main':
    solver = solve_system()
    # Get the amplitude reflection and transmission coefficients
    (rxCalculated, ryCalculated, rzCalculated) = (solver.rx, solver.ry, solver.rz)
    (txCalculated, tyCalculated, tzCalculated) = (solver.tx, solver.ty, solver.tz)

    # Get the diffraction efficiencies R and T and overall reflection and transmission coefficients R and T
    (R, T, RTot, TTot) = (solver.R, solver.T, solver.RTot, solver.TTot)
    print(RTot, TTot, RTot+TTot)
