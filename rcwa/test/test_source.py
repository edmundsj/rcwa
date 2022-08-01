"""
Module for testing the plane wave source class
"""
from rcwa import Layer, Source
from rcwa.shorthand import *
from numpy.testing import assert_allclose
import pytest
from copy import copy

@pytest.fixture
def source():
    wavelength = 0.02
    theta = 60 * deg
    phi = 30 * deg
    pTEM = 1 / sqrt(2) * complexArray([1, 1j])
    reflectionLayer = Layer(er=2, ur=1)
    src = Source(wavelength, theta, phi, pTEM, reflectionLayer)
    yield src

@pytest.mark.unit
def testKIncident(source):
    kIncidentActual = complexArray([1.0607, 0.61237, 0.70711])
    kIncidentCalculated = source.k_incident
    assert_allclose(kIncidentActual, kIncidentCalculated, atol=1e-5, rtol=1e-4)


@pytest.mark.unit
def test_equal_source(source):
    source_copy = copy(source)
    assert source == source_copy


@pytest.mark.unit
def test_equal_source_fail(source):
    with pytest.raises(ValueError):
        assert source == 2


@pytest.mark.unit
def test_print_source(source):
    assert len(str(source)) > 40


@pytest.mark.unit
def test_set_pTEM(source):
    source.pTEM = [2, 0]
    assert_allclose(source.pTEM , [1, 0])


@pytest.mark.unit
def test_set_phi(source):
    old_kincident = source.k_incident
    kinc_desired = np.array([np.sqrt(3)/2 * np.sqrt(2),  0, 1 / np.sqrt(2)])
    source.phi = 0
    assert old_kincident is not source.k_incident
    assert_allclose(source.k_incident, kinc_desired)