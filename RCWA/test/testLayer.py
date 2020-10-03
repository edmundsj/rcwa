import context

from crystal import *
from layer import *
from source import *
import unittest
from shorthandTest import *

class Test(unittest.TestCase):
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
