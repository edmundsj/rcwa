# Author: Jordan Edmunds, Ph.D. Student, UC Berkeley
# Contact: jordan.e@berkeley.edu
# Creation Date: 11/01/2019
#
import context
import RCWA

import numpy as np
import scipy as sp
import scipy.linalg
from matplotlib import pyplot as plt

lambda0 = 1.3
nSi = 3.5
nSiO2 = 1.45
tSi = lambda0/4/nSi
tSiO2 = lambda0/4/nSiO2
layer0 = RCWA.Layer(n=nSi)
layer1 = RCWA.Layer(n=nSiO2, L=tSiO2)
layer2 = RCWA.Layer(n=nSi, L=tSi)
layer3 = RCWA.Layer(n=nSiO2, L=tSiO2)
layer4 = RCWA.Layer(n=nSi, L=tSi)
layer5 = RCWA.Layer(n=nSiO2, L=tSiO2)
layer6 = RCWA.Layer(n=1)
stack = RCWA.LayerStack(layer0, layer1, layer2, layer3, layer4, layer5, layer6)
source = RCWA.Source(wavelength=0.4)

startWavelength = 0.4
stopWavelength = 1.5
stepWavelength = 0.002

print("Solving system...")
TMMSolver = RCWA.Solver(stack, source, (1, 1))
wavelengths = np.arange(startWavelength, stopWavelength + stepWavelength,
        stepWavelength)

TMMSolver.Solve(wavelengths=wavelengths)
#Plotter.plotEllipsometrySpectra(TMMSolver1.results)
RCWA.Plotter.plotRTEMSpectra(TMMSolver.results)
plt.show()
