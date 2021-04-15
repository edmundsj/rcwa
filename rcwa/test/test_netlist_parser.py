import context

import unittest
from shorthand import *
import numpy as np
import scipy as sp
import scipy.linalg
from matrices import *
from netlist_parser import *
from shorthandTest import *

# 1. The class NetlistParser parses a netlist and turns everything into a "Mask" or "Field" object.
# The masks and field are returned so that they are sorted in ascending order with
# respect to their coordinate on the optical axis.
class TestNetlistParser(unittest.TestCase):
    def setUp(self):
        self.netlist1arg = context.testLocation + '/netlists/test_netlist1arg.txt'
        self.netlist3arg =  context.testLocation + '/netlists/test_netlist3arg.txt'
        self.netlist4arg =  context.testLocation + '/netlists/test_netlist4arg.txt'
        self.netlist5arg =  context.testLocation + '/netlists/test_netlist5arg.txt'
        self.netlist6arg =  context.testLocation + '/netlists/test_netlist6arg.txt'
        self.netlist7arg =  context.testLocation + '/netlists/test_netlist7arg.txt'

    def testStripUnits(self):
        # First, test some simple real numbers
        parser = NetlistParser(self.netlist5arg)
        testString = "1um"
        valueActual = 1
        valueCalculated = parser.stripUnits(testString)
        assertEqual(valueActual, valueCalculated, errorMessage="testNetlistParser: testStripUnits: 1um")

        testString = "1.5nm"
        valueActual = 1.5e-3
        valueCalculated = parser.stripUnits(testString)
        assertEqual(valueActual, valueCalculated, errorMessage="testNetlistParser: testStripUnits: 1.5nm")

        testString = "1 + 4.0j"
        valueActual = 1 + 4.0j
        valueCalculated = parser.stripUnits(testString)
        assertEqual(valueActual, valueCalculated, errorMessage="testNetlistParser: testStripUnits: complex")

    def testProcessLines(self):
        actualLines = ['W0 1.4um 45deg 30deg 1 1j', 'LR 2', 'L0 1 1.5um', 'LT 3', '.REFRACTIVEINDEX']
        parser = NetlistParser(self.netlist5arg)
        parser.processLines()
        calculatedLines = parser.processedLines
        assertArrayEqual(actualLines, calculatedLines)

    def testExtractDirectives(self):
        refractiveIndexDirective = True
        parser = NetlistParser(self.netlist5arg)
        parser.processLines()
        parser.extractDirectives()

        actualLines = ['W0 1.4um 45deg 30deg 1 1j', 'LR 2', 'L0 1 1.5um', 'LT 3']
        calculatedLines = parser.processedLines
        assertArrayEqual(actualLines, calculatedLines,
                errorMessage="testNetlistParser: extractDirectives: lines")

        actualIndexDirective = True
        calculatedIndexDirective = parser.useRefractiveIndex
        assertEqual(actualIndexDirective, calculatedIndexDirective,
                errorMessage="testNetlistParser: extractDirectives: n")

    def testExtractLayers(self):
        parser = NetlistParser(self.netlist5arg)
        parser.processLines()
        parser.extractDirectives()
        parser.extractLayers()

        actualLines = ['W0 1.4um 45deg 30deg 1 1j']
        calculatedLines = parser.processedLines
        assertArrayEqual(actualLines, calculatedLines,
                errorMessage="testNetlistParser: extractLayers: lines")

        actualReflectionLayer = Layer(sq(2), 1)
        actualTransmissionLayer = Layer(sq(3), 1)
        actualInnerLayer = Layer(sq(2), 1, 1.5)
        actualLayerStack = LayerStack([actualReflectionLayer, actualInnerLayer, actualTransmissionLayer])
        calculatedLayerStack = parser.layerStack
        assertEqual(actualLayerStack, calculatedLayerStack, errorMessage='testNetlistParser: extractLayers: stacks unequal')

    def testExtractSources(self):
        parser = NetlistParser(self.netlist5arg)
        parser.processLines()
        parser.extractDirectives()
        parser.extractLayers()
        parser.extractSources()

        actualLines = []
        calculatedLines = parser.processedLines
        assertArrayEqual(actualLines, calculatedLines,
                errorMessage="testNetlistParser: extractSources: lines")

        actualReflectionLayer = Layer(sq(2), 1)
        actualSource = Source(1.4, 45*deg, 30*deg, [1, 1j], actualReflectionLayer)
        calculatedSource = parser.sources[0]
        assertAlmostEqual(actualSource, calculatedSource, errorMessage="testNetlistParser; extractSources: source")

        parser = NetlistParser(self.netlist1arg)
        parser.parseNetlist()
        actualSource = Source(1.4, 0*deg, 0*deg, [1, 0], actualReflectionLayer)
        calculatedSource = parser.sources[0]
        assertAlmostEqual(actualSource, calculatedSource, errorMessage="testNetlistParser; extractSources: source1")

        parser = NetlistParser(self.netlist3arg)
        parser.parseNetlist()
        actualSource = Source(1.4, 0*deg, 0*deg, [1, 0], actualReflectionLayer)
        calculatedSource = parser.sources[0]
        assertAlmostEqual(actualSource, calculatedSource, errorMessage="testNetlistParser; extractSources: source3")

        parser = NetlistParser(self.netlist4arg)
        parser.parseNetlist()
        actualSource = Source(1.4, 45*deg, 0*deg, [1, 1j], actualReflectionLayer)
        calculatedSource = parser.sources[0]
        assertAlmostEqual(actualSource, calculatedSource, errorMessage="testNetlistParser; extractSources: source4")

        parser = NetlistParser(self.netlist5arg)
        parser.parseNetlist()
        actualSource = Source(1.4, 45*deg, 30*deg, [1, 1j], actualReflectionLayer)
        calculatedSource = parser.sources[0]
        assertAlmostEqual(actualSource, calculatedSource, errorMessage="testNetlistParser; extractSources: source5")

        parser = NetlistParser(self.netlist6arg)
        parser.parseNetlist()
        actualSource = Source(1.4, 45*deg, 0*deg, [1, 1j], actualReflectionLayer)
        calculatedSource = parser.sources[0]
        assertAlmostEqual(actualSource, calculatedSource, errorMessage="testNetlistParser; extractSources: source6")

        parser = NetlistParser(self.netlist7arg)
        parser.parseNetlist()
        actualSource = Source(1.4, 45*deg, 30*deg, [1, 1j], actualReflectionLayer)
        calculatedSource = parser.sources[0]
        assertAlmostEqual(actualSource, calculatedSource, errorMessage="testNetlistParser; extractSources: source7")

    def testParseNetlist(self):
        parser = NetlistParser(self.netlist5arg)
        parser.parseNetlist()

        actualLines = []
        calculatedLines = parser.processedLines
        assertArrayEqual(actualLines, calculatedLines,
                errorMessage="testNetlistParser: extractSources: lines")

        actualReflectionLayer = Layer(sq(2), 1)
        actualSource = Source(1.4, 45*deg, 30*deg, [1, 1j], actualReflectionLayer)
        calculatedSource = parser.sources[0]
        assertAlmostEqual(actualSource, calculatedSource, errorMessage="testNetlistParser; parseNetlist: source")

        actualReflectionLayer = Layer(sq(2), 1)
        actualTransmissionLayer = Layer(sq(3), 1)
        actualInnerLayer = Layer(sq(2), 1, 1.5)
        actualLayerStack = LayerStack([actualReflectionLayer, actualInnerLayer, actualTransmissionLayer])
        calculatedLayerStack = parser.layerStack
        assertEqual(actualLayerStack, calculatedLayerStack, errorMessage='testNetlistParser: parseNetlist: stack')
