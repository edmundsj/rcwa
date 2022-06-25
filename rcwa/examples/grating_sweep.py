from rcwa import Source, Layer, LayerStack, Crystal, Solver, RectangularGrating
from rcwa.shorthand import complexArray
import numpy as np


def solve_system():
    reflection_layer = Layer(er=1.0, ur=1.0)
    transmission_layer = Layer(er=9.0, ur=1.0)

    wavelength = 0.5
    source = Source(wavelength=wavelength)

    N_harmonics = 11

    grating_layer = RectangularGrating(period=2, thickness=0.5, n=4, n_void=1, nx=500)
    layer_stack = LayerStack(grating_layer, incident_layer=reflection_layer, transmission_layer=transmission_layer)

    solver_1d = Solver(layer_stack, source, N_harmonics)
    solver_1d.solve((grating_layer, {'thickness': [0.3, 0.4, 0.5]}))

    return solver_1d

if __name__ == 'main':
    solver = solve_system()

    # Get the diffraction efficiencies R and T and overall reflection and transmission coefficients R and T
    (R, T, RTot, TTot) = (solver.R, solver.T, solver.RTot, solver.TTot)
    print(RTot, TTot, RTot+TTot)
