[![Build](https://github.com/edmundsj/rcwa/actions/workflows/build.yml/badge.svg)](https://github.com/edmundsj/rcwa/actions/workflows/build.yml) [![codecov](https://codecov.io/gh/edmundsj/rcwa/branch/master/graph/badge.svg?token=UDJ1TUESG3)](https://codecov.io/gh/edmundsj/rcwa) [![docs](https://github.com/edmundsj/rcwa/actions/workflows/build-docs.yml/badge.svg)](https://github.com/edmundsj/rcwa/actions/workflows/build-docs.yml) [![PyPI version](https://badge.fury.io/py/rcwa.svg)](https://badge.fury.io/py/rcwa) [![DOI](https://zenodo.org/badge/236611452.svg)](https://zenodo.org/badge/latestdoi/236611452)

What this package can do
===========================
- Calculate reflectance, transmittance, amplitude reflection coefficients, and scattering parameters from stacks of planar thin films
- Simulate reflectance, transmittance, diffraction efficiencies, scattering matrices from 1D diffraction gratings
- Simulate reflectance, transmittance, diffraction efficiencies, scattering matrices from 2D photonic crystals
- All these can be calculated fon incident plane waves with any polarization (linear, circular, elliptical), angle of incidence (phi and theta) and wavelength


Getting Started
================
Installation
--------------
The recommended way to install this software is with `pip`:

```
pip install rcwa
```

And that's it! 

Hello World Program
----------------------
To see a simple example, run:

```
python -m rcwa.examples.bragg_mirror
```

This should run an example with a 10-layer bragg mirror (also known as a [dielectric mirror](https://en.wikipedia.org/wiki/Dielectric_mirror)), which can have very high reflectance near its design wavelength, and output the reflectance as a function of wavelength, as seen below:

![Bragg Mirror Plot](/images/rcwa_example_plot.png)


Features
==========
- Implements 1D Transfer Matrix Method for homogenous layers
- Implements full rectangular 2D RCWA for periodic layers
- Huge material database for n/k values in optical range built-in based on [refractiveindex.info](https://refractiveindex.info/), including metals, plastics, glass, and ceramics
- Easy to use class-based syntax 
- Integrated parameter sweeps of any simulation parameter: geometry, materials, wavelength, angle of incidence, etc.
- Compute reflection and transmission spectra at arbitrary incidence and polarization
- Compute spectroscopic ellipsometry curves
- Compute reflected power, transmitted power, and S-parameters
- Large, fast-to-run test suite
- Extremely fast narrowband, rigorously correct simulations well suited for resonant devices
- Built-in convergence testing 

Example Uses
==============
- Compute reflected and transmitted power from a thin film stack
- Determine resonant frequency of a VCSEL
- Determine reflectance of a bragg mirror, on or off-axis
- Find diffraction efficiencies for a 1D or 2D diffraction grating
- Compute reflected power from a metallic mirror

Examples
============
All examples are in the `examples/` directory in your locally installed `rcwa` package, or in `rcwa/examples/` on this repository.

Reflection off Dispersive Materials
---------------------------------------

The below example demonstrates the reflection spectra you get reflecting off a bare surface of silicon, using the built-in materials database.

```
from rcwa import Material, Layer, LayerStack, Source, Solver, Plotter

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
startWavelength = 0.25
stopWavelength = 0.8
stepWavelength = 0.001

# Setup the source
source = Source(wavelength=startWavelength)

# Setup the materials and geometry
si = Material(name='Si')

# Setup the interface
reflectionLayer = Layer(n=1) # Free space
transmissionLayer = Layer(material=si)
stack = LayerStack(incident_layer=reflectionLayer, transmission_layer=transmissionLayer)

# Setup the solver
TMMSolver = Solver(stack, source, (1, 1))

# Setup and run the sweep
wavelengths = np.arange(startWavelength, stopWavelength + stepWavelength,
        stepWavelength)
results = TMMSolver.solve(wavelength=wavelengths)
results.plot(x='wavelength', y='RTot', show=True)
```
![Dispersive Si Plot](/images/si_dispersive.png)

Source Wavelength / Angle Sweeps
----------------------------------
```
import numpy as np
from rcwa import Material, Layer, LayerStack, Source, Solver, Plotter

# Setup the source
startWavelength = 0.25
stopWavelength = 0.8
stepWavelength = 0.02
wavelengths = np.arange(startWavelength, stopWavelength + stepWavelength,
        stepWavelength)
thetas = np.linspace(0, np.pi/4,10)

source = Source(wavelength=startWavelength)

thin_film = Layer(thickness=0.1, n=2)
substrate = Layer(n=4)
stack = LayerStack(thick_film, transmission_layer=substrate)

solver = Solver(stack, source)

results = solver.solve(wavelength=wavelengths, theta=thetas)
results.plot(x='wavelength', y='RTot', show=True)

```
![Reflectance vs Wavelength with varying angle](/images/wavelength_angle_sweep.png)

Geometry Sweeps
-------------------------------------------------------------------------------
Here, we set up a simulation with a rectangular grating on a substrate with a relative permittivity of 9, and a wavelength of 0.5 units (microns, meters, whatever you like!). This can be found in the `grating_sweep.py` example. In this example we are sweeping the thickness, but we could have swept the period or refractive index or permittivity.

```
from rcwa import Source, Layer, LayerStack, Crystal, Solver, RectangularGrating
import numpy as np
from matplotlib import pyplot as plt

reflection_layer = Layer(er=1.0, ur=1.0)
transmission_layer = Layer(er=9.0, ur=1.0)

wavelength = 0.5
source = Source(wavelength=wavelength)

N_harmonics = 11

grating_layer = RectangularGrating(period=2, thickness=0.5, n=4, n_void=1, nx=500)
layer_stack = LayerStack(grating_layer, incident_layer=reflection_layer, transmission_layer=transmission_layer)

solver_1d = Solver(layer_stack, source, N_harmonics)
results = solver_1d.solve((grating_layer, {'thickness': np.linspace(0.3, 0.5, 100)}))

results.plot(x='thickness', y='RTot', show=True)
```
![Reflectance vs Thickness](/images/reflectance_vs_thickness.png)

Documentation
================
This  project is documented on [Github Pages](https://edmundsj.github.io/rcwa/). For additional information, including downloading examples, you can view this project on [github](https://github.com/edmundsj/RCWA). 

Author: Jordan Edmunds, UC Irvine Alumnus, UC Berkeley Ph.D. Student

Date Started: 2020/01/05

What this package can't do
=========================================
- Simulate structures which are finite in x and y
- Simulate structures which are non-periodic or random

What this package could do with minor changes
==================================================
- Simulate responses to non-plane-wave incident electromagnetic fields
- Simulate modes inside gratings and photonic crystals and near-field patterns outside them
- Calculate diffraction orders from a 1D grating or 2D photonic crystal
- Use anisotropic materials
- Calculate band structures of photonic crystals or gratings

Frequently Asked Questions
=============================
Q: How do I tell the solver to use the Transfer Matrix Method or Rigorous Coupled Wave Analysis?

A: Don't worry, it will figure it out for you. If you define only unpatterned films, the solver will internally use the transfer matrix method. Otherwise, it will use RCWA.

Q: Will this work for tilted plane waves? What about elliptically polarized light?

A: Yes! RCWA and the TMM have no trouble with incident waves at arbitrary angles of incidence or polarization.

License
=========
This project is distributed under the [MIT license](https://mit-license.org/).

Dependencies
=============
Dependencies are comprehensively covered by the setup.py file, and the most
recent set of dependencies can be found there. Currently, this requires numpy,
scipy, pandas, matplotlib, and pyyaml. The documentation is built using Sphinx
and hosted on readthedocs.io.

Contributors
================
Contributors are welcome!  There's lots of interesting ways to expand this package, and lots of ideas on the "Issues" page. If you want to work on those, or something else, fork the repo and get started! :)

Acknowledgements / References
===============================
This work is based primarily on a set of lectures and associated course
material by Professor [Raymond Rumpf](http://emlab.utep.edu/team.htm)  at
the University of Texas, El Paso. 

[1] Rakić, Aleksandar D., Aleksandra B. Djurišić, Jovan M. Elazar, and Marian L. Majewski. "Optical properties of metallic films for vertical-cavity optoelectronic devices." Applied optics 37, no. 22 (1998): 5271-5283.

