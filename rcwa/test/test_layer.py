from rcwa import Crystal, Layer, LayerStack, Source, freeSpaceLayer
from rcwa.testing import *
import numpy as np
import pytest
from numpy.testing import assert_equal
from rcwa.shorthand import complexArray, complexIdentity
from matplotlib import pyplot as plt

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

def test_set_Kx():
    layer1 = Layer()
    layer2 = Layer()
    layer3 = Layer()
    incident_layer = Layer()
    transmission_layer = Layer()
    stack = LayerStack(layer1, layer2, layer3, incident_layer=incident_layer, transmission_layer=transmission_layer)
    Kx = np.array([[1,  2], [3, 4]])
    stack.Kx = Kx
    assert stack.Kx is  Kx
    assert layer1.Kx is Kx
    assert layer2.Kx is Kx
    assert layer3.Kx is Kx
    assert incident_layer.Kx is Kx
    assert transmission_layer.Kx is Kx
    assert stack.gapLayer.Kx is Kx

def test_set_Ky():
    layer1 = Layer()
    layer2 = Layer()
    layer3 = Layer()
    incident_layer = Layer()
    transmission_layer = Layer()
    stack = LayerStack(layer1, layer2, layer3, incident_layer=incident_layer, transmission_layer=transmission_layer)
    Ky = np.array([[1, 2], [3, 4]])
    stack.Ky = Ky
    assert stack.Ky is Ky
    assert layer1.Ky is Ky
    assert layer2.Ky is Ky
    assert layer3.Ky is Ky
    assert incident_layer.Ky is Ky
    assert transmission_layer.Ky is Ky
    assert stack.gapLayer.Ky is Ky


def test_all_layers():
    layer1 = Layer()
    i_layer = Layer()
    t_layer = Layer()
    stack = LayerStack(layer1, incident_layer=i_layer, transmission_layer=t_layer)
    all_layers = stack.all_layers
    assert len(all_layers) == 3
    assert all_layers[0] is i_layer
    assert all_layers[1] is layer1
    assert all_layers[2] is t_layer


def test_set_gap_layer():
    layer_i = Layer()
    layer_t = Layer()
    layer = Layer()
    stack = LayerStack(layer, incident_layer=layer_i, transmission_layer=layer_t)
    stack.Kx = 1.0006267892837541
    stack.Ky = 0.42474087247562836
    stack.set_gap_layer()
    WGap = complexIdentity(2)
    VGap = complexArray([
        [0 - 0.4250j, 0 - 1.1804j],
        [0 + 2.0013j, 0 + 0.4250j]])
    for layer in [layer_i, layer_t, layer]:
        assert_almost_equal(WGap, layer.Wg, absoluteTolerance=1e-4)
        assert_almost_equal(VGap, layer.Vg, absoluteTolerance=1e-4)


def  test_set_source():
    layer_i = Layer()
    layer_t = Layer()
    layer = Layer()
    source = Source()
    stack = LayerStack(layer, incident_layer=layer_i, transmission_layer=layer_t)
    stack.source = source
    assert stack.source is source
    assert stack.gapLayer.source is source
    for layer in [layer_i, layer_t, layer]:
        assert layer.source is source


def test_print_layer():
    layer = Layer()
    assert len(str(layer)) > 0


def test_print_stack():
    stack = LayerStack()
    assert len(str(stack)) > 0


def test_plot_newfigs():
    stack = LayerStack()
    fig, ax = stack.plot()
    assert fig is not None
    assert ax is not None


def test_plot_noax():
    stack = LayerStack()
    fig, ax = plt.subplots()
    new_fig, new_ax = stack.plot(fig=fig)
    assert new_fig is fig
    assert new_ax is not None
    assert new_ax is not ax


def test_plot_yesax():
    stack = LayerStack()
    fig, ax = plt.subplots()
    new_fig, new_ax = stack.plot(fig=fig, ax=ax)
    assert new_fig is fig
    assert new_ax is ax

