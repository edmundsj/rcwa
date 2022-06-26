import numpy as np
from rcwa import Material, Layer, LayerStack, Source, Solver, Plotter
import matplotlib.pyplot as plt

def solve_system():
    startWavelength = 0.25
    stopWavelength = 0.8
    stepWavelength = 0.02

    # Setup the source
    source = Source(wavelength=startWavelength)

    # Setup the Geometry
    thin_film = Layer(thickness=0.1, n=2)
    substrate = Layer(n=4)
    stack = LayerStack(thin_film, transmission_layer=substrate)

    # Setup the Solver
    solver = Solver(stack, source)

    # Setup and run the sweep
    wavelengths = np.arange(startWavelength, stopWavelength + stepWavelength,
            stepWavelength)
    thetas = np.linspace(0, np.pi/4,10)

    results = solver.solve(wavelength=wavelengths, theta=thetas)
    return results

if __name__ == '__main__':
    results = solve_system()
    angles, wavelengths, R = results['theta'], results['wavelength'], results['RTot']
    plt.plot(wavelengths, R)
    plt.show()
