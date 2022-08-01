from rcwa import Source, Layer, LayerStack, Crystal, Solver, RectangularGrating
from rcwa.shorthand import complexArray
import numpy as np


def solve_system():

    reflection_layer = Layer(er=1.0, ur=1.0)
    transmission_layer = Layer(er=9.0, ur=1.0)

    wavelength = 0.5
    deg = np.pi / 180
    k0 = 2*np.pi/wavelength
    theta = 60 * deg
    phi = 1*deg
    pTEM = 1/np.sqrt(2)*complexArray([1,1j])
    source = Source(wavelength=wavelength, theta=theta, phi=phi, pTEM=pTEM, layer=reflection_layer)

    crystal_thickness = 0.5

    N_harmonics = 11

    grating_layer = RectangularGrating(period=2, thickness=0.5, n=4, n_void=1, nx=500)
    layer_stack = LayerStack(grating_layer, incident_layer=reflection_layer, transmission_layer=transmission_layer)

    solver_1d = Solver(layer_stack, source, N_harmonics)
    results = solver_1d.solve()

    return results

if __name__ == '__main__':
    results = solve_system()

    # Get the amplitude reflection and transmission coefficients
    (rxCalculated, ryCalculated, rzCalculated) = (results['rx'], results['ry'], results['rz'])
    (txCalculated, tyCalculated, tzCalculated) = (results['tx'], results['ty'], results['tz'])

    # Get the diffraction efficiencies R and T and overall reflection and transmission coefficients R and T
    (R, T, RTot, TTot) = (results['R'], results['T'], results['RTot'], results['TTot'])
    print(RTot, TTot, RTot+TTot)
