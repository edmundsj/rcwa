[![Build](https://github.com/edmundsj/rcwa/actions/workflows/python-package-conda.yml/badge.svg)](https://github.com/edmundsj/rcwa/actions/workflows/python-package-conda.yml) [![codecov](https://codecov.io/gh/edmundsj/rcwa/branch/master/graph/badge.svg?token=UDJ1TUESG3)](https://codecov.io/gh/edmundsj/rcwa) [![docs](https://github.com/edmundsj/rcwa/actions/workflows/build-docs.yml/badge.svg)](https://github.com/edmundsj/rcwa/actions/workflows/build-docs.yml) [![PyPI version](https://badge.fury.io/py/rcwa.svg)](https://badge.fury.io/py/rcwa) [![DOI](https://zenodo.org/badge/236611452.svg)](https://zenodo.org/badge/latestdoi/236611452)


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
To run a simple example, run:

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
- Compute reflection and transmission spectra at arbitrary incidence and polarization
- Compute spectroscopic ellipsometry curves
- Exactly solves Maxwell's Equations for arbitrary layer stacks of any thickness
- Compute reflected power, transmitted power, and S-parameters
- Easy to use class-based syntax 
- Large, fast-to-run test suite
- Extremely fast narrowband, rigorously correct simulations well suited for resonant devices

Example Uses
==============
- Compute reflected and transmitted power from a thin film stack
- Determine resonant frequency of a VCSEL
- Determine reflectance of a bragg mirror, on or off-axis
- Find diffraction efficiencies for a 1D or 2D diffraction grating
- Compute reflected power from a metallic mirror

Documentation
================
This  project is documented on [Github Pages](https://edmundsj.github.io/rcwa/). For additional information, including downloading examples, you can view this project on [github](https://github.com/edmundsj/RCWA). 

Author: Jordan Edmunds, UC Irvine Alumnus, UC Berkeley Ph.D. Student

Date Started: 2020/01/05

License
=========
This project is distributed under the [MIT license](https://mit-license.org/).

Dependencies
=============
Dependencies are comprehensively covered by the setup.py file, and the most
recent set of dependencies can be found there. Currently, this requires numpy,
scipy, pandas, matplotlib, and pyyaml. The documentation is built using Sphinx
and hosted on readthedocs.io.

Acknowledgements / References
===============================
This work is based primarily on a set of lectures and associated course
material by Professor [Raymond Rumpf](http://emlab.utep.edu/team.htm)  at
the University of Texas, El Paso. 

[1] Rakić, Aleksandar D., Aleksandra B. Djurišić, Jovan M. Elazar, and Marian L. Majewski. "Optical properties of metallic films for vertical-cavity optoelectronic devices." Applied optics 37, no. 22 (1998): 5271-5283.

