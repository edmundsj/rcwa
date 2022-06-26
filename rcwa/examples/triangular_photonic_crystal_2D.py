from rcwa import Source, Layer, LayerStack, Crystal, Solver, example_dir
from rcwa.shorthand import complexArray
import numpy as np
import os


def solve_system():
    data_path = os.path.join(example_dir, 'triangleData.csv')
    devicePermittivityCellData = np.transpose(np.loadtxt(data_path, delimiter=','))
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

    deviceCrystal = Crystal(t1, t2, er=devicePermittivityCellData, ur=devicePermeabilityCellData)
    layer1 = Layer(crystal=deviceCrystal, thickness=crystalThickness)
    layerStack = LayerStack(layer1, incident_layer=reflectionLayer, transmission_layer=transmissionLayer)

    rcwa_solver = Solver(layerStack, source, numberHarmonics)
    results = rcwa_solver.solve()
    return results


if __name__ == '__main__':
    results = solve_system()
    # Get the amplitude reflection and transmission coefficients
    (rxCalculated, ryCalculated, rzCalculated) = (results['rx'], results['ry'], results['rz'])
    (txCalculated, tyCalculated, tzCalculated) = (results['tx'], results['ty'], results['tz'])

    # Get the diffraction efficiencies R and T and overall reflection and transmission coefficients R and T
    (R, T, RTot, TTot) = (results['R'], results['T'], results['RTot'], results['TTot'])
    print(RTot, TTot, RTot + TTot)
