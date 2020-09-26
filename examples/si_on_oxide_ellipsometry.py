# Author: Jordan Edmunds, Ph.D. Student, UC Berkeley
# Contact: jordan.e@berkeley.edu
# Creation Date: 11/01/2019
#
# TODO:

import numpy as np
import scipy as sp
import scipy.linalg
import sys
from RCWA.source.source import Source
from RCWA.source.layer import Layer, LayerStack
from RCWA.source.solver import RCWASolver
from RCWA.netlist.netlist_parser import *
from RCWA.source.plotter import Plotter

import matplotlib.pyplot as plt

layer0 = Layer(n=1)
layer1 = Layer(n=3.5, L=0.092)
layer2 = Layer(n=1.45, L=0.226)
layer3 = Layer(n=3.5, L=0.092)
layer4 = Layer(n=1.45, L=0.226)
layer5 = Layer(n=3.5)

layerStack = LayerStack(layer0, layer1, layer2, layer3, layer4, layer5)
source = Source(wavelength=0.5, theta=np.radians(0))
startWavelength = 0.4
stopWavelength = 1.7
stepWavelength = 0.001

TMMSolver = RCWASolver(layerStack=layerStack, source=source, numberHarmonics=(1,1))

print("Solving system...")
wavelengths = np.arange(startWavelength, stopWavelength + stepWavelength,
        stepWavelength)

TMMSolver.Solve(wavelengths=wavelengths)
#Plotter.plotEllipsometrySpectra(TMMSolver.results)
Plotter.plotReflectionSpectra(TMMSolver.results)
plt.show()
print("Done!")
