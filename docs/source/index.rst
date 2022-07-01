.. Rigorous Coupled Wave Analysis documentation master file, created by
   sphinx-quickstart on Mon Sep 28 12:56:28 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Documentation for Rigorous Coupled Wave Analysis
==========================================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   examples
   api
   material
   mathematics

Getting Started
------------------
Installation
______________
The recommended way to install this software is with `pip`:

.. code-block::

    pip install rcwa

And that's it!

Hello World Program
----------------------
To run a simple example, run:

.. code-block::

    python -m rcwa.examples.bragg_mirror


This should run an example with a 10-layer bragg mirror (also known as a [dielectric mirror](https://en.wikipedia.org/wiki/Dielectric_mirror)), which can have very high reflectance near its design wavelength, and output the reflectance and transmission as a function of wavelength.

2D Photonic Crystal Example
------------------------------
The examples folder also contains a 2D photonic crystal example called `triangular_photonic_crystal.py`. Running it will print the reflection and transmission coefficients, but diving into the example you can extract any quantities of interest: amplitude reflection coefficients, diffraction efficiencies for each harmonic, etc.

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

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
