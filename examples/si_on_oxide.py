# Author: Jordan Edmunds, Ph.D. Student, UC Berkeley
# Contact: jordan.e@berkeley.edu
# Creation Date: 11/01/2019
#
import sys

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

layer0 = Layer(n=3.5)
layer1 = Layer(n=1.45, L=0.94)
layer2 = Layer(n=3.5, L=0.03)
layer3 = Layer(n=1)
stack = LayerStack(layer0, layer1, layer2, layer3)
source = Source(wavelength=0.4)

startWavelength = 0.4
stopWavelength = 1
stepWavelength = 0.02

print("Solving system...")
TMMSolver1 = RCWASolver(stack, source, (1, 1))
wavelengths = np.arange(startWavelength, stopWavelength + stepWavelength,
        stepWavelength)

TMMSolver1.Solve(wavelengths=wavelengths)
#Plotter.plotEllipsometrySpectra(TMMSolver1.results)
Plotter.plotRTEMSpectra(TMMSolver1.results)
plt.show()
