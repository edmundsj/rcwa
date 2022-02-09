# Author: Jordan Edmunds, Ph.D. Student, UC Berkeley
# Contact: jordan.e@berkeley.edu
# Creation Date: 11/01/2019
#
import context
from rcwa import Material, Layer, LayerStack, Source, Solver, Plotter

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

designWavelength = 1.3
startWavelength = 0.25
stopWavelength = 0.85
stepWavelength = 0.001

source = Source(wavelength=designWavelength)
siO2 = Material(material_name='SiO2')

reflectionLayer = Layer(n=1) # Free space
transmissionLayer = Layer(material=siO2)
stack = LayerStack(reflectionLayer, transmissionLayer)

print("Solving system...")
TMMSolver = Solver(stack, source, (1, 1))
wavelengths = np.arange(startWavelength, stopWavelength + stepWavelength,
        stepWavelength)

TMMSolver.Solve(wavelengths=wavelengths)
#Plotter.plotEllipsometrySpectra(TMMSolver.results)
Plotter.plotRTSpectra(TMMSolver.results)
plt.show()
