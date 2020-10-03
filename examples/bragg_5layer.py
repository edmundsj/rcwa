# Author: Jordan Edmunds, Ph.D. Student, UC Berkeley
# Contact: jordan.e@berkeley.edu
# Creation Date: 11/01/2019
#
import sys
import RCWA

import numpy as np
import scipy as sp
import scipy.linalg
sys.path.append('../core/')
#sys.path.append('../')
from layer import Layer, LayerStack
from source import Source
from plotter import Plotter
from matplotlib import pyplot as plt
#from RCWA import RCWASolver # this is the syntax I want
from solver import RCWASolver # This is the syntax I currently have

lambda0 = 1.3
nSi = 3.5
nSiO2 = 1.45
tSi = lambda0/4/nSi
tSiO2 = lambda0/4/nSiO2
layer0 = Layer(n=nSi)
layer1 = Layer(n=nSiO2, L=tSiO2)
layer2 = Layer(n=nSi, L=tSi)
layer3 = Layer(n=nSiO2, L=tSiO2)
layer4 = Layer(n=nSi, L=tSi)
layer5 = Layer(n=nSiO2, L=tSiO2)
layer6 = Layer(n=1)
stack = RCWA.LayerStack(layer0, layer1, layer2, layer3, layer4, layer5, layer6)
source = RCWA.Source(wavelength=0.4)

startWavelength = 0.4
stopWavelength = 1.5
stepWavelength = 0.002

print("Solving system...")
TMMSolver1 = RCWASolver(stack, source, (1, 1))
wavelengths = np.arange(startWavelength, stopWavelength + stepWavelength,
        stepWavelength)

TMMSolver1.Solve(wavelengths=wavelengths)
#Plotter.plotEllipsometrySpectra(TMMSolver1.results)
Plotter.plotRTEMSpectra(TMMSolver1.results)
plt.show()
