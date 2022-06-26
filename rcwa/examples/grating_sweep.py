from rcwa import Source, Layer, LayerStack, Crystal, Solver, RectangularGrating
import numpy as np
from matplotlib import pyplot as plt


def solve_system():
    reflection_layer = Layer(er=1.0, ur=1.0)
    transmission_layer = Layer(er=9.0, ur=1.0)

    wavelength = 0.5
    source = Source(wavelength=wavelength)

    N_harmonics = 11

    grating_layer = RectangularGrating(period=2, thickness=0.5, n=4, n_void=1, nx=500)
    layer_stack = LayerStack(grating_layer, incident_layer=reflection_layer, transmission_layer=transmission_layer)

    solver_1d = Solver(layer_stack, source, N_harmonics)
    results = solver_1d.solve((grating_layer, {'thickness': np.linspace(0.3, 0.5, 100)}))

    return results


if __name__ == '__main__':
    results = solve_system()
    plt.plot(results['thickness'], results['RTot'])
    plt.show()