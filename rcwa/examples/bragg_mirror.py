# Author: Jordan Edmunds, Ph.D. Student, UC Berkeley
# Contact: jordan.e@berkeley.edu
# Creation Date: 11/01/2019
#
import context
import rcwa

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

reflectionLayer = rcwa.Layer(n=1)
transmissionLayer = rcwa.Layer(n=n1)
layer0 = rcwa.Layer(n=n1, L=t1)
layer1 = rcwa.Layer(n=n2, L=t2)
layer2 = rcwa.Layer(n=n1, L=t1)
layer3 = rcwa.Layer(n=n2, L=t2)
layer4 = rcwa.Layer(n=n1, L=t1)
layer5 = rcwa.Layer(n=n2, L=t2)
layer6 = rcwa.Layer(n=n1, L=t1)
layer7 = rcwa.Layer(n=n2, L=t2)
layer8 = rcwa.Layer(n=n1, L=t1)
layer9 = rcwa.Layer(n=n2, L=t2)
layer10 = rcwa.Layer(n=n1, L=t1)
stack = rcwa.LayerStack(reflectionLayer,
       layer0, layer1, layer2, layer3, layer4, layer5, layer6, layer7, layer8, layer9, layer10,
        transmissionLayer)
source = rcwa.Source(wavelength=designWavelength)

print("Solving system...")
TMMSolver = rcwa.Solver(stack, source, (1, 1))
wavelengths = np.arange(startWavelength, stopWavelength + stepWavelength,
        stepWavelength)

TMMSolver.Solve(wavelengths=wavelengths)
#Plotter.plotEllipsometrySpectra(TMMSolver1.results)
rcwa.Plotter.plotRTSpectra(TMMSolver.results)
plt.show()
