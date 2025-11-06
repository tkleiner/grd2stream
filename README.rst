grd2stream
==========

**grd2stream** is a command-line utility that performs 4th-order Runge–Kutta streamline
integration from gridded 2D velocity fields. The program outputs streamline coordinates
that can directly be used with the Generic Mapping Tools (`GMT <https://www.generic-mapping-tools.org/>`_)
or QGIS for scientific visualization and mapping.

Usage
-----

Basic usage example:

.. code-block:: bash

    # bash, GMT4
    echo "0 0" | grd2stream vx.grd vy.grd | psxy -m -R -J ...

Newer versions of grd2stream can be linked to the
`GMT-API <https://docs.generic-mapping-tools.org/latest/devdocs/api.html>`_ for file IO. The old GMT5-API is no longer supported.

If GMT6 is installed with GDAL support, all GDAL file formats (e.g., NetCDF, GTiFF, etc.) can be read.
If you are not familiar with the GMT-style of accessing the different file formats, please check
the `GMT reference <https://docs.generic-mapping-tools.org/latest/reference/features.html#write-grids-images>`_.

.. code-block:: bash

    # bash, GMT6
    echo "0 0" | grd2stream velocity.nc?vx velocity.nc?vy | gmt plot ...

Note, in In tcsh, the question mark (``?``) is a wildcard that matches exactly one character.
If you want to use a literal '?', you must escape it so the shell doesn’t treat it as a wildcard.

.. code-block:: tcsh

    # tcsh, GMT6
    echo "0 0" | grd2stream velocity.nc\?vx velocity.nc\?vy | gmt plot ...


Example
-------

.. figure:: /docs/_static/map_FRIS_B13_line_source_small.png
   :scale: 50 %
   :alt: map of Filchner-Ronne Ice Shelf (FRIS) including ice core location B13

   Map of Filchner-Ronne Ice Shelf (FRIS) including ice core location B13. 
   Several streamlines generated from point or line sources.


Installation
------------

**Installing required libraries on macOS**

These instructions assume that you use `MacPorts <https://www.macports.org/>`_.
Note, you can have gmt4 and gmt6 installed at the same time without interference.

.. code-block:: console

    sudo port install gdal +netcdf+hdf5
    sudo port install gmt6 +gdal

**Installing required libraries using Conda**

.. code-block:: console

    # on albedo0|1
    module purge
    module load conda
    # setup your environment called GMT6 as an example
    conda create -y -n GMT6 gmt=6* gdal hdf5 netcdf4

**Install from Release Tarball**

Before you start to compile grd2stream, make sure gmt-config is in your ``PATH``.

.. code-block:: console

  tar xvfz grd2stream-X.X.X.tar.gz
  cd grd2stream-X.X.X
  ./configure
  make && make install

By default, grd2stream is installed in your ``$HOME`` folder.
Make sure ``$HOME/bin`` is in your ``PATH``.

You can specify a custom installation prefix:

.. code-block:: console

  ./configure --prefix=/usr/local
  make
  sudo make install

Alternatively install grd2stream into your GMT6 conda environment

.. code-block:: console

  conda activate GMT6
  ./configure --prefix=$CONDA_PREFIX
  make && make install

If GMT is not available on your machine, grd2stream can still be used with
with the the classic GMT grid file format (netCDF2).

If gmt-config or nc-config can not be found by configure, you could specify the paths

.. code-block:: console

  ./configure --with-gmt-config=/opt/local/lib/gmt6/bin/gmt-config --with-nc-config=/opt/local/bin/nc-config


**Install from Git repository**

This requires **GNU Autotools** (autoconf, automake, aclocal, libtool, and make) to be installed

.. code-block:: console

  git clone https://github.com/tkleiner/grd2stream.git
  cd grd2stream
  ./bootsrap.sh
  ./configure
  make && make install

Packaging (Maintainer only)
---------------------------

Build a release tarball with:

.. code-block:: bash

  git clone https://github.com/tkleiner/grd2stream.git
  cd grd2stream
  ./bootsrap.sh
  ./configure
  make && make dist

To remove all auto-generated file run ``make total-clean``.
It deletes files that require bootstrapping to rebuild.


Fourth-order Runge–Kutta (RK4) Streamline Integration
-----------------------------------------------------

A streamline in a 2D steady velocity field :math:`\mathbf{v}(x, y) = (v_x(x, y), v_y(x, y))`
is a curve :math:`\mathbf{s}(t) = (x(t), y(t))` that satisfies the ordinary differential
equation (ODE):

.. math::

    \frac{d\mathbf{s}}{dt} = \mathbf{v}(\mathbf{s}(t)) =
    (v_x(x(t), y(t)), v_y(x(t), y(t)))

with initial condition:

.. math::

    \mathbf{s}(t_0) = \mathbf{s}_0 = (x_0, y_0)

Here, :math:`t` is a parametric variable along the streamline.

To numerically approximate :math:`\mathbf{s}(t)`, the classical 4th-order Runge–Kutta
scheme with fixed step size :math:`h > 0` is used.

Given :math:`\mathbf{s}_i = (x_i, y_i)` at step :math:`i`, the next point
:math:`\mathbf{s}_{i+1}` is computed as:

.. math::

    \mathbf{k}_1 &= \mathbf{v}\left(\mathbf{s}_i\right), \\
    \mathbf{k}_2 &= \mathbf{v}\left(\mathbf{s}_i + \frac{h}{2}\mathbf{k}_1\right), \\
    \mathbf{k}_3 &= \mathbf{v}\left(\mathbf{s}_i + \frac{h}{2}\mathbf{k}_2\right), \\
    \mathbf{k}_4 &= \mathbf{v}\left(\mathbf{s}_i + h\,\mathbf{k}_3\right), \\
    \mathbf{s}_{i+1} &= \mathbf{s}_i + \frac{h}{6}\left(\mathbf{k}_1 + 2\mathbf{k}_2 + 2\mathbf{k}_3 + \mathbf{k}_4\right) = \mathbf{s}_i + \mathbf{\delta}_i

Component-wise, this can be written explicitly as:

.. math::

    k_{1x} &= v_x(x_i, y_i), \quad & k_{1y} &= v_y(x_i, y_i), \\
    k_{2x} &= v_x(x_i + \frac{h}{2}k_{1x}, y_i + \frac{h}{2}k_{1y}), \quad &
    k_{2y} &= v_y(x_i + \frac{h}{2}k_{1x}, y_i + \frac{h}{2}k_{1y}), \\
    k_{3x} &= v_x(x_i + \frac{h}{2}k_{2x}, y_i + \frac{h}{2}k_{2y}), \quad &
    k_{3y} &= v_y(x_i + \frac{h}{2}k_{2x}, y_i + \frac{h}{2}k_{2y}), \\
    k_{4x} &= v_x(x_i + h k_{3x}, y_n + h k_{3y}), \quad &
    k_{4y} &= v_y(x_i + h k_{3x}, y_n + h k_{3y}),

with

.. math::

    x_{i+1} &= x_i + \frac{h}{6}(k_{1x} + 2 k_{2x} + 2 k_{3x} + k_{4x}) = x_i + dx_i, \\
    y_{i+1} &= y_i + \frac{h}{6}(k_{1y} + 2 k_{2y} + 2 k_{3y} + k_{4y}) = y_i + dy_i .

The integration of a streamline continues until one of the following
conditions is met:

- The number of integration steps exceeds a preset maximum :math:`N_{\max}`.
- The streamline point leaves the computational domain.
- The velocity magnitude :math:`|\mathbf{v}(x, y)|` becomes smaller than a threshold,
  indicating a stagnation point.

Acknowledgments
---------------

* Generic Mapping Tools (`GMT <https://www.generic-mapping-tools.org/>`_)
* `GDAL <https://gdal.org/>`_ by the Open Source Geospatial Foundation (`OSGeo <https://www.osgeo.org/>`_)

Contribute
----------

- Issue Tracker: https://github.com/tkleiner/grd2stream/issues
