# Author: Jordan Edmunds, Ph.D. Student, UC Berkeley
# Contact: jordan.e@berkeley.edu
# Creation Date: 11/01/2019
#
import context
from rcwa import Layer, LayerStack, Source, Solver, Plotter

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

reflectionLayer = Layer(n=1)
transmissionLayer = Layer(n=n1)
layer0 = Layer(n=n1, L=t1)
layer1 = Layer(n=n2, L=t2)
layer2 = Layer(n=n1, L=t1)
layer3 = Layer(n=n2, L=t2)
layer4 = Layer(n=n1, L=t1)
layer5 = Layer(n=n2, L=t2)
layer6 = Layer(n=n1, L=t1)
layer7 = Layer(n=n2, L=t2)
layer8 = Layer(n=n1, L=t1)
layer9 = Layer(n=n2, L=t2)
layer10 = Layer(n=n1, L=t1)
stack = LayerStack(reflectionLayer,
       layer0, layer1, layer2, layer3, layer4, layer5, layer6, layer7, layer8, layer9, layer10,
        transmissionLayer)
source = Source(wavelength=designWavelength)

print("Solving system...")
TMMSolver = Solver(stack, source, (1, 1))
wavelengths = np.arange(startWavelength, stopWavelength + stepWavelength,
        stepWavelength)

TMMSolver.Solve(wavelengths=wavelengths)
#Plotter.plotEllipsometrySpectra(TMMSolver1.results)
Plotter.plotRTSpectra(TMMSolver.results)
plt.show()
