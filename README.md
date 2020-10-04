[![Build Status](https://travis-ci.com/edmundsj/RCWA.svg?branch=master)](https://travis-ci.com/edmundsj/RCWA) [![Coverage Status](https://coveralls.io/repos/github/edmundsj/RCWA/badge.svg?branch=master)](https://coveralls.io/github/edmundsj/RCWA?branch=master) [![Documentation Status](https://readthedocs.org/projects/rcwa/badge/?version=latest)](https://rcwa.readthedocs.io/en/latest/?badge=latest)[![PyPI version](https://badge.fury.io/py/rcwa.svg)](https://badge.fury.io/py/rcwa)

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

This should run an example with a 10-layer bragg mirror (also known as a [dielectric mirror](https://en.wikipedia.org/wiki/Dielectric_mirror)), which can have very high reflectance near its design wavelength, and output the reflectance as a function of wavelength.

Features
==========
- Implements 1D Transfer Matrix Method for homogenous layers
- Implements full rectangular 2D RCWA for periodic layers
- Arbitrary incident wave polarization (circular, linear, elliptical)
- Arbitrary incident wave angle of incidence
- Exactly solves Maxwell's Equations for arbitrary layer stacks of any thickness
- Compute reflected power, transmitted power, and S-parameters
- Easy to use class-based syntax 
- Large, fast-to-run test suite

Documentation
================
This  project is documented on [Read The Docs](https://rcwa.readthedocs.io/en/latest/). For additional information, including downloading examples, you can view this project on [github](https://github.com/edmundsj/RCWA). 

Author: Jordan Edmunds, UC Irvine Alumnus, UC Berkeley Ph.D. Student

Date Started: 2020/01/05

License
=========
This project is distributed under the [MIT license](https://mit-license.org/).

