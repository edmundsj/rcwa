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
        stopWavelength = 0.85
        stepWavelength = 0.001

        source = Source(wavelength=startWavelength)
        si = Material(name='Si')
        data = pd.DataFrame({'Wavelength (um):': si.wavelengths, 'er': si._er_dispersive, 'n': si._n_dispersive})
        print(data)

        reflectionLayer = Layer(n=1) # Free space
        transmissionLayer = Layer(material=si)
        stack = LayerStack(reflectionLayer, transmissionLayer)

        print("Solving system...")
        TMMSolver = Solver(stack, source, (1, 1))
        wavelengths = np.arange(startWavelength, stopWavelength + stepWavelength,
                stepWavelength)
        TMMSolver.solve(wavelength=wavelengths)
        Plotter.plotRTSpectra(TMMSolver.results)
        return TMMSolver

if __name__ == 'main':
        solver = solve_system()
        plt.show()
