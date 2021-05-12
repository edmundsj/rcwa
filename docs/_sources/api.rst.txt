.. Rigorous Coupled Wave Analysis documentation master file, created by
   sphinx-quickstart on Mon Sep 28 12:56:28 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

rcwa API
==========================================================

There are four main user-facing classes: :code:`Solver`, :code:`Layer`, :code:`LayerStack`, :code:`Source`, :code:`Material`, :code:`Crystal`, and :code:`Plotter`. The :code:`Layer` class is used to set up the individual simulation layers and define their properties, including thickness, material properties, and crystal structures (if applicable). The :code:`Material` class can be passed into the :code:`Layer` class to define a layer whose material properties are dependent on refractive index. The :code:`Crystal` class is used to define layers which are periodic in x and y, and is also passed into the :code:`Layer` class. Once all individual layers are constructed, they are used to create a :code:`LayerStack`, which contains information about which region is considered the "incident" region and which is the "transmission" region (both are semi-infinite).

The :code:`Source` class is used to define the excitation source - the wavelength, polarization, and incident angle. 

Once the user has created a :code:`Source` and :code:`LayerStack` class, these are passed into the :code:`Solver` class, which then runs the simulation, and makes the results available as a dictionary :code:`Solver.results`.

.. autoclass:: rcwa.Solver
    :members: Solve

.. autoclass:: rcwa.Layer

.. autoclass:: rcwa.LayerStack

.. autoclass:: rcwa.Source

.. autoclass:: rcwa.Crystal

.. autoclass:: rcwa.Material

.. autoclass:: rcwa.Plotter
    :members: plotRTSpectra

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
