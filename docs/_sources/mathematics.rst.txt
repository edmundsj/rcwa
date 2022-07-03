.. Rigorous Coupled Wave Analysis documentation master file, created by
   sphinx-quickstart on Mon Sep 28 12:56:28 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Mathematics of RCWA
==========================================================

This section is for developers who want to understand the implementation of this package, for example, to implement their own manipulations on top of it or extend functionality. It follows the formulation laid out by `Rumpf <https://empossible.net/wp-content/uploads/2019/08/Lecture-7a-RCWA-Formulation.pdf>`_. This is intended as a reference for developers with a fairly advanced electromagnetics and mathematics background, and assumes a knowledge of electromagnetic modes, vectors, matrices, matrix multiplication, and vectors and matrices composed of other vectors and matrices.

Definitions and Conventions
-----------------------------

   | :math:`N_{x}`: number of harmonics along x-direction
   | :math:`N_y`: number of harmonics along y-direction.
   | :math:`N_{tot} =  N_x * N_y`: Total number of harmonics
   | :math:`n_x`: Used to index the x-harmonics. Ranges from negative to positive.
   | :math:`n_y`: Used to index the y-harmonics. Ranges from negative to positive.
   | :math:`E`: Electric Field
   | :math:`H`: Magnetic Field
   | :math:`s`: Electtric field (Fourier) coefficient
   | :math:`u`: Magnetic field (Fourier) coefficient
   | :math:`c`: Mode coefficient
   | :math:`W`: Electric field coefficient mode matrix
   | :math:`V`: Magnetic field coefficient mode matrix

Layer "0" is the incident layer. Layer "1" is the layer closest to the incident region.

Fields and Field Coefficients
----------------------------------
The formulation of rigorous coupled wave analysis involves decomposing the electric and magnetic fields into their plane wave components. Rather than working with the electric field :math:`E_x` directly,  we work instead with the field coefficients :math:`s_x`.

.. math::
   E_x  (x, y, z) = \sum_{n_x = -\infty}^{\infty} \sum_{n_y= -\infty}^{\infty} s_x(n_x, n_y, z) * e^{-j \left(k_x(n_x, n_y)*x + k_y(n_x, n_y)*y \right)}

Where

.. math::

   k_x(n_x, n_y) = k_{x, incident} - n_x * T_{1, x} - n_y * T_{2, x} \\
   k_y(n_x, n_y) = k_{y, incident} - n_x * T_{1, y} - n_y * T_{2, y}

where :math:`T_1` and :math:`T_2` are the lattice parameters of the photonic crystal under study. If a 1D grating is used, :math:`T_2 = 0`. For thin films, :math:`T_1` is also zero, and only one harmonic needs to be represented - the incident harmonic.

Note that the field coefficients are, in general, a function of z.

Field Coefficients
------------------------
The electric (:math:`s`) and magnetic (:math:`u`) field coefficients are internally represented by vectors which have both x- and y-components, and these are different for each. Taking just electric field coefficient vector in the i-th layer:

.. math::

   s_{i}(z)  = \begin{bmatrix}
           s_{i, x}(z) \\
            s_{i, y}(z)
   \end{bmatrix}

The length of the :math:`s_{i,x}` and :math:`s_{i,y}` vectors have length equal to :math:`N_{tot}`, so the total length of the :math:`s_{i}` vector is :math:`N_{tot}*2`. These vectors are indexed :math:`(n_x, n_y)`. For a :math:`N_x=3` and :math:`N_y=3` setup, the vector looks like this:

.. math::

   s_{i, x}  = \begin{bmatrix}
   s_{i, x}(1, 1) \\
   s_{i, x}(1, 2) \\
   s_{i, x}(1, 3) \\
   s_{i, x}(2, 1) \\
   s_{i, x}(2, 2) \\
   s_{i, x}(2, 3) \\
   s_{i, x}(3, 1) \\
   s_{i, x}(3, 2) \\
   s_{i, x}(3, 3) \\
   \end{bmatrix}

The x-component of the zero-order harmonic corresponds to the field coefficient :math:`s_{i,x}(2, 2)`.

These field coefficients are, in general, functions of the z-coordinate. This makes them  awkward to work with directly, which is why we typically work with *mode coefficients* instead.

Mode Coefficients
-------------------------------------------------------

Internally, this package typically works with *mode coefficients*, referred to using the symbol :math:`c`.  Forward-propagating mode coefficients inside the :math:`i^{th}` layer are represented with the symbol :math:`c_{i}^+` and backward-propagating mode coefficients inside the :math:`i^{th}` layer are represented with the symbol :math:`c_{i}^-`. Each mode has a single value of :math:`k_z`, the propagation constant in the :math:`z`-direction. You could say this is what *defines* a mode. For planar films, the modes are always plane waves. For 1D- and 2D-photonic crystals, they require several plane waves (harmonics) to represent.

Just as with field coefficients, both the x- and y-components of the mode coefficients must be kept internally. As with the field coefficients, both forward- and backward-propagating mode coefficients are represented with the x-components followed by the y-components.

.. math::

   c_{i}^+  = \begin{bmatrix}
           c_{i, x}^+ \\
            c_{i, y}^+
   \end{bmatrix}

:math:`c_{i,x}^+` and :math:`c_{i,y}^+` are vectors with length equal to the total number of harmonics (i.e. for a 2D photonic crystal if :math:`N_x = 3` and :math:`N_y=3`, then :math:`N_{tot} = N_x * N_y = 9`, and the length of the total vector :math:`c_{i}^+` is 18.

Scattering Matrices couple Mode Coefficients
------------------------------------------------
The scattering matrix for the :math:`i^{th}` layer relates the forward- and backward propagating mode coefficients between the :math:`i^{th}` and :math:`i+1^{th}` layer in the following way:

.. math::

   \begin{bmatrix}
           c_{i}^- \\
            c_{i+1}^+
   \end{bmatrix} = \begin{bmatrix}
      S_{i,11} & S_{i, 12}\\
      S_{i, 21} & S_{i, 22}
   \end{bmatrix} \begin{bmatrix}
      c_{i}^+ \\
      c_{i+1}^-
   \end{bmatrix}

This can be rearranged to solve for the :math:`i+1^{th}` mode coefficients given the :math:`i^{th}` mode coefficients.

Electric Field Coefficients
______________________________

The mode coefficients in the :math:`i^{th}` layer can be related to directly to the x- and y-components of the electric fields (:math:`s_x` and :math:`s_y`) using the electric mode matrix, represented with the symbol :math:`W`. This matrix takes the forward- and backward-traveling mode coefficients and converts them into the electric field coefficients. The values of :math:`k_z` for each mode are also needed, which are assembled diagonally along the :math:`\lambda` matrix.

.. math::

   \begin{bmatrix}
           s_{i,x} \\
            s_{i,y}
   \end{bmatrix} = W . e^{- \lambda z} . c_{i}^+ + W . e^{- \lambda z} . c_{i}^-

Magnetic Field Coefficients
___________________________________________

Similarly, the mode coefficients in the :math:`i^{th}` layer can be related to the x- and y-components of the magnetic field coefficients :math:`u_x` and :math:`u_y` using the magnetic mode matrix, represented with the symbol :math:`V`, along with, as before, the :math:`\lambda` matrix.

.. math::

   \begin{bmatrix}
           u_{i,x} \\
            u_{i,y}
   \end{bmatrix} = - V . e^{- \lambda z} . c_{i}^+ + V . e^{- \lambda z} . c_{i}^-

Finding Mode coefficients inside an arbitrary layer
-------------------------------------------------------
Once the scattering matrices for each layer :math:`S_i` are known, the incident mode coefficients :math:`s_0` are known and the global scattering matrix :math:`S_{global}` is known, the mode coefficients in an arbitrary layer can be calculated.

First, the mode coefficients in the incident region must be found. To do this, you can find the m. :math:`c_0^+` using the :math:`W_{incident}` matrix and :math:`s_{incident}` vector (which contains the x- and y-components of the zero-order harmonic), and :math:`c_0^-` can be calculated from the global scattering matrix. Note that :math:`s_{incident}` is not the same as :math:`s_0`, which also contains the reflected field coefficients.

.. math::
   c_0^+ = W_{incident}^{-1} s_{incident}
   c_0^- = S_{global, 11} c_0^+

Then, by applying the formula below as many times as is required, mode coefficients within the desired layer can be found:

.. math::

   c_{i+1}^{-} = S_{i, 12}^{-1} c_i^{-} - S_{i, 12}^{-1} S_{i,11} c_i{+} \\
   c_{i+1}^{+} = S_{i, 21} c_{i}^{+1} + S_{i,22} c_{i+1}{-}

Find E/H Coefficients inside an arbitrary Layer
------------------------------------------------------
Once the mode coefficients have been found, the electric and magnetic field coefficients can be found as described previously. Note that at this point, the field coefficients will still be a function of the z coordinate.

Finding the electric and magnetic fields inside an arbitrary layer
-----------------------------------------------------------------------

Electric Fields
_________________
Once the electric field coefficients are known, the electric fields can be calculated using the formula above for a certain value of :math:`x` and :math:`y`, as a function of :math:`z`.

First, the diagonal k-matrices :math:`K_x` and :math:`K_y` can be converted into vectors :math:`k_x` and :math:`k_y`. Then, the complex exponential can be evaluated element-wise directly:

.. math::

   e_{vec} = e^{-j \left(k_x * x + k_y * y \right) }

This is itself a vector with the same length as :math:`s_x`. The dot product (NOT inner product, there is no complex conjugation) of this vector and the exponential vector :math:`e_{vec}` finally yields the x-component of the electric field. This may be repeated for the y-component as well, by dotting with the :math:`s_y` vector.

.. math::

   E_x(x, y, z) = e_{vec} \cdot s_{x} \\
