# Author: Jordan Edmunds, Ph.D. Student, UC Berkeley
# Contact: jordan.e@berkeley.edu
# Creation Date: 11/01/2019
#
from rcwa import Material, Layer, LayerStack, Source, Solver, Plotter
import numpy as np
from matplotlib import pyplot as plt

def solve_system():
        startWavelength = 0.25
        stopWavelength = 0.850
        stepWavelength = 0.001
        incident_angle = np.radians(75)

        source = Source(wavelength=startWavelength, theta=incident_angle)
        si = Material('Si')

        reflectionLayer = Layer(n=1) # Free space
        transmissionLayer = Layer(material=si)
        stack = LayerStack(incident_layer=reflectionLayer, transmission_layer=transmissionLayer)

        TMMSolver = Solver(stack, source)
        wavelengths = np.arange(startWavelength, stopWavelength + stepWavelength,
                stepWavelength)

        TMMSolver.solve(wavelength=wavelengths)

        tan_psi_predicted = TMMSolver.results['tanPsi']
        cos_delta_predicted = TMMSolver.results['cosDelta']

        fig, ax = plt.subplots()
        ax.plot(wavelengths, tan_psi_predicted)
        ax.plot(wavelengths, cos_delta_predicted)
        ax.legend([r'$Tan(\Psi)$', r'$Cos(\Delta)$'])

        return fig, ax


if __name__ == '__main__':
        fig, ax = solve_system()
        plt.show()
