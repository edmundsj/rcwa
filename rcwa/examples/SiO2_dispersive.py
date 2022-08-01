# Author: Jordan Edmunds, Ph.D. Student, UC Berkeley
# Contact: jordan.e@berkeley.edu
# Creation Date: 11/01/2019
#
from rcwa import Material, Layer, LayerStack, Source, Solver, Plotter
import numpy as np
from matplotlib import pyplot as plt

def solve_system():
        startWavelength = 0.25
        stopWavelength = 0.85
        stepWavelength = 0.001

        source = Source(wavelength=startWavelength)
        siO2 = Material('SiO2')

        reflectionLayer = Layer(n=1) # Free space
        transmissionLayer = Layer(material=siO2)
        stack = LayerStack(incident_layer=reflectionLayer, transmission_layer=transmissionLayer)

        TMMSolver = Solver(stack, source, (1, 1))
        wavelengths = np.arange(startWavelength, stopWavelength + stepWavelength,
                stepWavelength)

        results = TMMSolver.solve(wavelength=wavelengths)
        return results


if __name__ == '__main__':
        results = solve_system()
        fig, ax = results.plot(x='wavelength', y=['RTot', 'TTot', 'conservation'])
        plt.show()
