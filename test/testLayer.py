import sys
sys.path.append('core');
sys.path.append('test')

from crystal import *
from layer import *
import unittest
from shorthandTest import *

class Test(unittest.TestCase):
    def testExtractCrystal(self):
        testCrystal = Crystal()
        testLayer = Layer(crystal=testCrystal)
        testStack = LayerStack(freeSpaceLayer, testLayer, freeSpaceLayer)

        internalLayerActual = 0
        internalLayerCalculated = testStack.extractCrystalLayer()
        assertEqual(internalLayerActual, internalLayerCalculated)

        testStack = LayerStack(freeSpaceLayer, freeSpaceLayer, testLayer, freeSpaceLayer)

        internalLayerActual = 1
        internalLayerCalculated = testStack.extractCrystalLayer()
        assertEqual(internalLayerActual, internalLayerCalculated)
