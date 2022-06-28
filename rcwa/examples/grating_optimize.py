from autograd import grad

from rcwa import Source, Layer, LayerStack, Crystal, Solver, RectangularGrating
import numpy as np
from matplotlib import pyplot as plt

def solve_inner(thickness):
    reflection_layer = Layer(er=1.0, ur=1.0)
    transmission_layer = Layer(er=9.0, ur=1.0)

    wavelength = 0.5
    source = Source(wavelength=wavelength)

    N_harmonics = 11

    grating_layer = RectangularGrating(period=2, thickness=thickness, n=4, n_void=1, nx=500)
    layer_stack = LayerStack(grating_layer, incident_layer=reflection_layer, transmission_layer=transmission_layer)

    solver_1d = Solver(layer_stack, source, N_harmonics)
    results = solver_1d.solve()
    return results

def loss_func(thickness):
    results = solve_inner(thickness)
    return results['RTot']


def solve_system():
    loss_grad = grad(loss_func)
    return loss_grad

if __name__ == '__main__':
    my_grad = solve_system()
    print(my_grad(0.5))
    #plt.plot(results['thickness'], results['RTot'])
    plt.show()