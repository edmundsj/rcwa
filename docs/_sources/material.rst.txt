.. Rigorous Coupled Wave Analysis documentation master file, created by
   sphinx-quickstart on Mon Sep 28 12:56:28 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Materials in rcwa
==========================================================

.. autoclass:: rcwa.Material

The `Material` class handles materials whose material properties (permittivity, permeability, or refractive index) change as a function of wavelength. There are three primary ways to describe materials:

**1. Constant numerical value**

   This is the simplest model for a material - a constant permittivity / permeability or refractive index, passed in with the `er`, `ur` and `n` arguments.

**2. Custom function of wavelength**

   In addition to being constants, `er`, `ur`, and `n` may be user-specified single-argument functions, which take wavelength as an argument. Units of wavelength in this case do not matter and returns a (potentially complex) scalar.

**3. Tabulated data**

    If a value is passed to the :code:`filename` argument, then the :code:`Material` will load in that file. Csv formats are supported by default, with the first axis being wavelength and the second axis being the complex-valued refractive index. You may instead choose to separate the real- and imaginary parts into two separate columns. This works as well. An example tabulated data file is shown below.

**4. Using a materials database - tabulated or dispersion formula**

    The refractiveindex.info database is supported by default, and contains many user-specified materials (i.e. :code:`Si`, :code:`SiO2`, :code:`Ti`, :code:`Pt`. This database provides both tabulated data and dispersion formulas (see `dispersion formulas <https://refractiveindex.info/database/doc/Dispersion%20formulas.pdf>`_), depending on the material used. CAUTION: if using this database, wavelength units must be specified in micrometers. 

Interpolation and Extrapolation
-------------------------------------
If using tabulated data, you may use any wavelength between the minimum and maximum values in the table, and the permittivity will be linearly interpolated if the desired wavelength falls between two tabulated points. If the desired wavelength falls outside the minimum or maximum, :code:`rcwa` will attempt to extrapolate using the slope near the boundary and give a warning.

Imaginary Sign Convention
----------------------------------------------------------------------
Lossy materials are represented with a positive imaginary refractive index (i.e. :code:`2 + 0.1j`. Materials with gain are represented by a negative imaginary refractive index (i.e. :code:`2 - 0.1j`).

Example Tabulated Data Files
----------------------------------------------------------------------

Below is what an example file for tabulated n/k data should look like, let's call it :code:`custom_material.csv`:

.. code-block::

    wavelength, n, k
    0.76, 2, 0.1
    0.77, 2.1, 0.2
    0.78, 2.1, 0.4

And a valid alternative :code:`custom_material2.csv`:

.. code-block::

    wavelength, n
    0.76, 2 + 0.1j
    0.77, 2.1 + 0.2j
    0.78, 2.1 + 0.4j

Note: the files must have a header.

Material Examples
---------------------------------

.. code-block::

    from rcwa import Material

    # Use a constant index or permittivity
    my_material = Material(n=5 + 0.1j)
    my_material_2 = Material(er=4, ur=1.1)

    # Use a custom function
    def n_func(wavelength):
        return 1 + wavelength / 1.5

    my_dispersive_material = Material(n=n_func)

    # Use built-in databases
    Si = Material('Si')
    SiO2 = Material('SiO2')

    # Use a custom file
    custom_material = Material(filename='custom_filename.csv')

