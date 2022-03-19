# Author: Jordan Edmunds, Ph.D. Student, UC Berkeley
# Contact: jordan.e@berkeley.edu
# Creation Date: 11/01/2019
#
import context
from rcwa import Material, Layer, LayerStack, Source, Solver, Plotter

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

startWavelength = 0.25
stopWavelength = 0.850
stepWavelength = 0.001
incident_angle = np.radians(75)

source = Source(wavelength=startWavelength, theta=incident_angle)
si = Material('Si')
breakpoint()

reflectionLayer = Layer(n=1) # Free space
transmissionLayer = Layer(material=si)
stack = LayerStack(reflectionLayer, transmissionLayer)

print("Solving system...")
TMMSolver = Solver(stack, source, (1, 1))
wavelengths = np.arange(startWavelength, stopWavelength + stepWavelength,
        stepWavelength)

TMMSolver.Solve(wavelengths=wavelengths)

tan_psi_predicted = np.array([result['tanPsi'] for result in TMMSolver.results])
cos_delta_predicted = np.array([result['cosDelta'] for result in TMMSolver.results])

fig, ax = plt.subplots()
ax.plot(wavelengths, tan_psi_predicted)
ax.plot(wavelengths, cos_delta_predicted)
ax.legend([r'$Tan(\Psi)$', r'$Cos(\Delta)$'])
plt.show()
