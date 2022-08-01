# Author: Jordan Edmunds, Ph.D. Student, UC Berkeley
# Contact: jordan.e@berkeley.edu
# Creation Date: 11/01/2019
#
from rcwa import Material, Layer, LayerStack, Source, Solver, Plotter

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

def solve_system():
        startWavelength = 0.25
        stopWavelength = 0.8
        stepWavelength = 0.001

        # Setup the source
        source = Source(wavelength=startWavelength)

        # Setup the materials and geometry
        si = Material(name='Si')

        # Setup the interface
        reflectionLayer = Layer(n=1) # Free space
        transmissionLayer = Layer(material=si)
        stack = LayerStack(incident_layer=reflectionLayer, transmission_layer=transmissionLayer)

        # Setup the solver
        TMMSolver = Solver(stack, source, (1, 1))

        # Setup and run the sweep
        wavelengths = np.arange(startWavelength, stopWavelength + stepWavelength,
                stepWavelength)
        results = TMMSolver.solve(wavelength=wavelengths)
        return results


if __name__ == '__main__':
        results = solve_system()
        results.plot(x='wavelength', y=['RTot', 'TTot', 'conservation'])
        plt.show()
