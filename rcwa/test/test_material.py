import pytest
import os
import numpy as np
from rcwa import Material, Source
from rcwa.testing import *
from rcwa.shorthand import *

@pytest.fixture
def silicon():
    source = Source(wavelength=1)
    silicon = Material(database_path='main/Si/Schinke.yml',source=source)
    yield silicon

def testExtractMaterialDatabase():
    fake_material = Material()
    materials_to_check = ['Pt', 'Si', 'Ag', 'Ti', 'Au', 'SiO2']
    assert all(i in fake_material.database.materials.keys() for i in materials_to_check) == True

def testLoadFromDatabase():
    wavelength = 0.1879 # Should load from Johnson nk table
    source = Source(wavelength=wavelength)
    Ag = Material('Ag', source=source)
    n_desired = 1.07+1.212j
    n_observed = Ag.n
    assert_almost_equal(n_desired, n_observed, absoluteTolerance=1e-3)

    source.wavelength=0.2262
    n_desired = 1.26+1.344j
    n_observed = Ag.n
    assert_almost_equal(n_desired, n_observed, absoluteTolerance=1e-3)

def testLoadFromFile():
    source = Source(wavelength=1.25)
    test_dir = os.path.dirname(os.path.abspath(__file__))
    filename = os.path.join(test_dir, 'data', 'ito.csv')
    desired_wavelengths = [1200, 1250, 1300]
    desired_n = 1.5 + 0.7j
    ito = Material(filename=filename, source=source)
    actual_n = ito.n
    assert actual_n == desired_n

def testnk(silicon):
    # Case 1: wavelength is exactly identical to one we have in database
    nDesired = 1.737 + 3.9932j
    silicon.source.wavelength = 0.26
    nCalculated = silicon.n
    assert_almost_equal(nCalculated, nDesired, errorMessage="material: testnk: n1", absoluteTolerance=1e-3)


    # Case 3: wavelength is larger than one we have in database - extrapolate linearly
    silicon.source.wavelength = 1.47
    nDesired = 3.485 + 1.09e-13j
    with pytest.warns(Warning) as warning:
        nCalculated = silicon.n

    assert_almost_equal(nCalculated, nDesired, errorMessage="material: testnk: n3", absoluteTolerance=1e-3)

    # Case 4: wavelength is smaller than one we have in database - extrapolate in the opposite direction
    silicon.source.wavelength = 0.23
    nDesired = 1.437 + 2.7803j
    with pytest.warns(Warning) as warning:
        nCalculated = silicon.n
    assert_almost_equal(nCalculated, nDesired, errorMessage="material: testnk: n4", absoluteTolerance=1e-3)

def testnkInterpolate(silicon):
    # Case 2: wavelength located between two points in the database - interpolate 
    nDesired = 4.2620000000000005 + 0.0461865j
    silicon.source.wavelength = 0.505
    nCalculated = silicon.n
    assert_almost_equal(nCalculated, nDesired, errorMessage="material: testnk: n21", absoluteTolerance=1e-3)

    # Case 2 repeated: wavelength located between two points in the database - interpolate 
    nDesired = 4.2485 + 0.0450088J
    silicon.source.wavelength = 0.5075
    nCalculated = silicon.n
    assert_almost_equal(nCalculated, nDesired, errorMessage="material: testnk: n22", absoluteTolerance=1e-3)

    # Case 2 software is having issues with: wavelength of 0.495 to 0.496. First 0.495
    nDesired = 4.319 + 0.051176j
    silicon.source.wavelength = 0.495
    nCalculated = silicon.n
    assert_almost_equal(nCalculated, nDesired, errorMessage="material: testnk: n23", absoluteTolerance=1e-3)

    # Now 0.496
    nDesired = 4.313 + 0.0506492j
    silicon.source.wavelength = 0.496
    nCalculated = silicon.n
    assert_almost_equal(nCalculated, nDesired, errorMessage="material: testnk: n24", absoluteTolerance=1e-3)

""" Tests for discontinuities in the Aspnes data, which I haven't found in Schinke for Si"""
@pytest.mark.unit
def testAvoidDiscontinuities():
    source = Source(wavelength = 0.495)
    silicon = Material(database_path='main/Si/Aspnes.yml', source=source)

    n_desired = 4.325778947368422 + 0.07380526315789473j
    n_observed = silicon.n
    assert_almost_equal(n_observed, n_desired, errorMessage="material: testnk: n25", absoluteTolerance=1e-6)

    source.wavelength = 0.496
    n_desired = 4.319492753623189 + 0.07293719806763285j
    n_observed = silicon.n
    assert_almost_equal(n_observed, n_desired, errorMessage="material: testnk: n26", absoluteTolerance=1e-6)


@pytest.mark.unit
def testEr(silicon):
    # Case 1: wavelength is exactly identical to one we have in database
    erDesired = sq(1.737 + 3.9932j)
    silicon.source.wavelength = 0.26
    erCalculated = silicon.er
    assert_almost_equal(erCalculated, erDesired, errorMessage="material: testnk: er1")

    # Case 2: wavelength interpolated between two points in the database
    silicon.source.wavelength = 0.264
    erDesired = -14.557277+15.787005j
    erCalculated = silicon.er
    assert_almost_equal(erCalculated, erDesired, errorMessage="material: testnk: er2", absoluteTolerance=1e-5)

    # Case 3: wavelength is larger than one we have in database - extrapolate linearly
    silicon.source.wavelength = 1.47
    erDesired = sq(3.487 + 1.09e-13j) - 0.01395
    with pytest.warns(Warning) as warning:
        erCalculated = silicon.er
    assert_almost_equal(erCalculated, erDesired, absoluteTolerance=1e-5, errorMessage="material: testnk: er3")

    # Case 4: wavelength is smaller than one we have in database
    silicon.source.wavelength = 0.23
    erDesired = -4.744348+7.505422j
    with pytest.warns(Warning) as warning:
        erCalculated = silicon.er
    assert_almost_equal(erCalculated, erDesired, absoluteTolerance = 1e-5, errorMessage="material: testnk: er4")


@pytest.mark.unit
def testUr(silicon):
    # Case 1: wavelength is exactly identical to one we have in database
    urDesired = 1
    silicon.wavelength = 0.27
    urCalculated = silicon.ur
    assert_almost_equal(urCalculated, urDesired, errorMessage="material: testnk: ur1")

    # Case 2: wavelength is nearly identical to one we have in database
    wavelength = 0.264
    urDesired = 1
    urCalculated = silicon.ur
    assert_almost_equal(urCalculated, urDesired, errorMessage="material: testnk: ur2")

    # Case 3: wavelength is larger than one we have in database - extrapolate linearly
    wavelength = 1.47
    urDesired = 1
    urCalculated = silicon.ur
    assert_almost_equal(urCalculated, urDesired, absoluteTolerance=1e-5, errorMessage="material: testnk: ur3")

    # Case 4: wavelength is smaller than one we have in database
    wavelength = 0.23
    urDesired = 1
    urCalculated = silicon.ur
    assert_almost_equal(urCalculated, urDesired, absoluteTolerance = 1e-5, errorMessage="material: testnk: ur4")


@pytest.mark.unit
def test_extract_dispersion_formula_2():
    src = Source(wavelength=0.5)
    SiO2 = Material(database_path='main/SiO2/Ghosh-e.yml', source=src)
    n_desired = 1.5580
    n_actual = SiO2.n
    assert_almost_equal(n_actual, n_desired, absoluteTolerance=1e-5)


@pytest.mark.unit
def test_extract_dispersion_formula_1():
    src = Source(wavelength=0.5)
    SiO2 = Material(database_path='main/SiO2/Radhakrishnan-o.yml', source=src)
    n_desired = 1.548755
    n_actual = SiO2.n
    assert_almost_equal(n_actual, n_desired, absoluteTolerance=1e-5)


@pytest.mark.unit
def test_dispersive_func():
    src = Source(wavelength=0.1)
    mat = Material(er=lambda x: x, ur=lambda x: 2*x, source=src)
    assert mat.er == 0.1
    assert mat.ur == 0.2

@pytest.mark.unit
def test_dispersive_func_er():
    src = Source(wavelength=0.1)
    mat = Material(er=lambda x: x, ur=2, source=src)
    assert mat.er == 0.1
    assert mat.ur == 2


@pytest.mark.unit
def test_dispersive_func_ur():
    src = Source(wavelength=0.1)
    mat = Material(er=3, ur=lambda x: x, source=src)
    assert mat.er == 3
    assert mat.ur == 0.1