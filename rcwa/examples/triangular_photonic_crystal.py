import unittest
from rcwa import Source, Layer, LayerStack, Crystal, Plotter, complexArray, Solver
import numpy as np

devicePermittivityCellData = np.transpose(np.loadtxt('triangleData.csv', delimiter=','))
devicePermeabilityCellData = 1 + 0 * devicePermittivityCellData

reflectionLayer = Layer(er=2.0, ur=1.0)
transmissionLayer = Layer(er=9.0, ur=1.0)

wavelength = 2
deg = np.pi / 180
k0 = 2*np.pi/wavelength
theta = 60 * deg
phi = 30*deg
pTEM = 1/np.sqrt(2)*complexArray([1,1j])
source = Source(wavelength=wavelength, theta=theta, phi=phi, pTEM=pTEM, layer=reflectionLayer)
t1, t2 = complexArray([1.75, 0, 0]), complexArray([0, 1.5, 0])

crystalThickness = 0.5

numberHarmonics = (3, 3)

deviceCrystal = Crystal(devicePermittivityCellData, devicePermeabilityCellData, t1, t2)
layer1 = Layer(crystal=deviceCrystal, L=crystalThickness, numberHarmonics=numberHarmonics)
layerStack = LayerStack(reflectionLayer, layer1, transmissionLayer)

solver = Solver(layerStack, source, numberHarmonics)
solver.Solve()

# Get the amplitude reflection and transmission coefficients
(rxCalculated, ryCalculated, rzCalculated) = (solver.rx, solver.ry, solver.rz)
(txCalculated, tyCalculated, tzCalculated) = (solver.tx, solver.ty, solver.tz)

# Get the diffraction efficiencies R and T and overall reflection and transmission coefficients R and T
(R, T, RTot, TTot) = (solver.R, solver.T, solver.RTot, solver.TTot)
print(RTot, TTot, RTot+TTot)
