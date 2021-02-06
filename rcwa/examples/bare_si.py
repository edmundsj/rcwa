# Author: Jordan Edmunds, Ph.D. Student, UC Berkeley
# Contact: jordan.e@berkeley.edu
# Creation Date: 11/01/2019
#
import context
from rcwa import Material, Layer, LayerStack, Source, Solver

import numpy as np
from matplotlib import pyplot as plt

designWavelength = 1.3
startWavelength = 0.6
stopWavelength = 2.2
stepWavelength = 0.005

n1 = 3.5 # refractive index of layer 1 (Si)
n2 = 1.45 # refractive index of layer 2 (SiO2)
t1 = designWavelength/4/n1
t2 = designWavelength/4/n2
si = Material('Si')

reflectionLayer = Layer(n=1)
transmissionLayer = Layer(material=si)
stack = LayerStack(reflectionLayer, transmissionLayer)
source = Source(wavelength=designWavelength)

print("Solving system...")
TMMSolver = Solver(stack, source, (1, 1))
wavelengths = np.arange(startWavelength, stopWavelength + stepWavelength,
        stepWavelength)

TMMSolver.Solve(wavelengths=wavelengths)
#Plotter.plotEllipsometrySpectra(TMMSolver1.results)
rcwa.Plotter.plotRTSpectra(TMMSolver.results)
plt.show()
