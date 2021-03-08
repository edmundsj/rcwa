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
si = Material('Pt')
bto = Material('Si')

reflectionLayer = Layer(n=1) # Free space
internal_layer = Layer(material=bto, L=0.1)
transmissionLayer = Layer(material=si)
stack = LayerStack(reflectionLayer, internal_layer, transmissionLayer)

print("Solving system...")
TMMSolver = Solver(stack, source, (1, 1))
wavelengths = np.arange(startWavelength, stopWavelength + stepWavelength,
        stepWavelength)

TMMSolver.Solve(wavelengths=wavelengths)

tan_psi_predicted = np.array([result['tanPsi'] for result in TMMSolver.results])
cos_delta_predicted = np.array([result['cosDelta'] for result in TMMSolver.results])

measured_data = pd.read_csv(context.nkLocation + '/Pt-ellipsometry-data.csv', delimiter='\t')
predicted_data = pd.DataFrame({'Wavelength (nm)': wavelengths * 1000, 'Tan(Psi)': tan_psi_predicted, 'Cos(Del)': -cos_delta_predicted})

fig, ax = plt.subplots()
measured_data.plot(x='Wavelength (nm)', y='Tan(Psi)',ax=ax, fig=fig)
predicted_data.plot(x='Wavelength (nm)', y='Tan(Psi)', ax=ax, fig=fig)
ax.legend([r'$Tan(\Psi) Measured$', r'$Tan(\Psi) Predicted$'])
plt.show()

fig, ax = plt.subplots()
measured_data.plot(x='Wavelength (nm)', y='Cos(Del)',ax=ax, fig=fig)
predicted_data.plot(x='Wavelength (nm)', y='Cos(Del)', ax=ax, fig=fig)
ax.legend([r'$Cos(\Delta) Measured$', r'$Cos(\Delta) Predicted$'])
plt.show()
