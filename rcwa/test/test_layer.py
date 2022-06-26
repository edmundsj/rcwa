from rcwa import Crystal, Layer, LayerStack, Source, freeSpaceLayer
from rcwa.testing import *
import numpy as np
import pytest
from numpy.testing import assert_equal

@pytest.fixture
def crystal_1D():
    permittivity_data = np.array([1, 1, 1, 1, 3, 3, 3, 3])
    permeability_data = 1 + 0 * permittivity_data
    t1 = np.array([1, 0])
    crystal = Crystal(permittivity_data, permeability_data, t1)
    return crystal

def testExtractCrystal():
    t1 = np.array([1, 0, 0])
    t2 = np.array([0, 1, 0])
    testCrystal = Crystal(t1, t2)
    testLayer = Layer(crystal=testCrystal)
    testStack = LayerStack(freeSpaceLayer, testLayer, freeSpaceLayer)

    actualCrystal = testCrystal
    calculatedCrystal = testStack.crystal
    assert_equal(testCrystal, actualCrystal)

    testStack = LayerStack(freeSpaceLayer, freeSpaceLayer, testLayer, freeSpaceLayer)

    calculatedCrystal = testStack.crystal
    assert_equal(calculatedCrystal, actualCrystal)

def test_empty_layer():
    empty_layer = Layer()
    desired_defaults = {'n': 1, 'er': 1, 'ur': 1, 'source': None}
    for key, val in desired_defaults.items():
        assert hasattr(empty_layer, key)
    for key, desired_val in desired_defaults.items():
        actual_val = getattr(empty_layer, key)
        assert actual_val == desired_val

def test_empty_stack():
    empty_stack = LayerStack()
    desired_defaults = {'internal_layers': [], 'incident_layer': freeSpaceLayer, 'transmission_layer': freeSpaceLayer}
    for key, val in desired_defaults.items():
        assert hasattr(empty_stack, key)
    for key, desired_val in desired_defaults.items():
        actual_val = getattr(empty_stack, key)
        assert actual_val == desired_val