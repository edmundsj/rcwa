import unittest

from RCWA.source.crystal import *
from RCWA.source.layer import *
from RCWA.source.source import *
from RCWA.source.material import Material
from RCWA.test.shorthandTest import *

# TODO - NEED TO WRITE TEST FOR SETWAVELENGTH FUNCTION.
# TODO - BREAK OUT TESTS FROM TESTMATRICES ETC. INTO THEIR OWN FILES

class LayerTest(unittest.TestCase):
    def testExtractCrystal(self):
        t1 = np.array([1,0,0])
        t2 = np.array([0,1,0])
        testCrystal = Crystal(t1, t2)
        testLayer = Layer(crystal=testCrystal)
        testStack = LayerStack(freeSpaceLayer, testLayer, freeSpaceLayer)

        internalLayerActual = 0
        internalLayerCalculated = testStack.extractCrystalLayer()
        assertEqual(internalLayerActual, internalLayerCalculated)

        testStack = LayerStack(freeSpaceLayer, freeSpaceLayer, testLayer, freeSpaceLayer)

        internalLayerActual = 1
        internalLayerCalculated = testStack.extractCrystalLayer()
        assertEqual(internalLayerActual, internalLayerCalculated)

    def testSetWavelength(self):
        materialLocation = 'nkData/Si_Schinke.csv'
        silicon = Material(materialLocation)
        testLayer = Layer(material=silicon)
        testLayer.setWavelength(0.25)
        desiredIndex = 1.637 + 3.59j
        desiredEr = sq(desiredIndex)
        desiredUr = 1
        calculatedIndex = testLayer.n
        calculatedEr = testLayer.er
        calculatedUr = testLayer.ur
        assertAlmostEqual(calculatedIndex, desiredIndex, errorMessage="Layer: setWavelength: n")
        assertAlmostEqual(calculatedEr, desiredEr, errorMessage="Layer: setWavelength: er")
        assertAlmostEqual(calculatedUr, desiredUr, errorMessage="Layer: setWavelength: ur")

