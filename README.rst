Implements a Runge-Kutta stream line integration based on gridded velocity
fields in native GMT grid format (vx.grd, vy.grd).

.. code-block:: bash

    echo "0 0" | grd2stream vx.grd vy.grd | psxy -m -R -J ...

Newer versions of `grd2stream` can be linked to the
`GMT-API <https://docs.generic-mapping-tools.org/latest/devdocs/api.html>`_ for file IO. The old GMT5-API is no longer supported.

If GMT6 is installed with GDAL support, all GDAL file formats (e.g., NetCDF, GTiFF, etc.) can be read.
If you are not familiar with the GMT-style of excessing the different file formats, please check
the `GMT reference <https://docs.generic-mapping-tools.org/6.5/reference/features.html#write-grids-images>`_.

.. code-block:: bash

    echo "0 0" | grd2stream velocity.nc?vx velocity.nc?vy | gmt plot ...


Example
-------

.. figure:: /docs/_static/map_FRIS_B13_line_source_small.png
   :scale: 50 %
   :alt: map of Filchner-Ronne Ice Shelf (FRIS) including ice core location B13

   Map of Filchner-Ronne Ice Shelf (FRIS) including ice core location B13. 
   Several stream lines as point or line sources from grd2stream and GMT mapping tools. 


Installation
------------

**Installing required libraries on macOS**

These instructions assume that you use `MacPorts <https://www.macports.org/>`_.
Note, you can have gmt4 and gmt6 installed at the same time without interference.

.. code-block:: bash

    sudo port install gdal +netcdf+hdf5
    sudo port install gmt6 +gdal

**Installing required libraries using Conda**

.. code-block:: bash

    # on albedo0|1
    module purge
    module load conda
    # setup your environment called GMT6 as an example
    conda create -y -n GMT6 gmt=6* gdal hdf5 netcdf4

**Install grd2stream version X.X.X**

.. code-block:: bash

  tar xvfz grd2stream-X.X.X.tar.gz
  cd grd2stream-X.X.X
  ./configure --prefix=$HOME --enable-gmt-api
  make && make install

Alternatively install `grd2stream` into your GMT6 conda environment

.. code-block:: bash

  conda activate GMT6
  tar xvfz grd2stream-X.X.X.tar.gz
  cd grd2stream-X.X.X
  ./configure ./configure --prefix=$CONDA_PREFIX --enable-gmt-api
  make && make install


Packaging (Maintainer only)
---------------------------

Build grd2stream-X.X.X.tar.gz e.g. by

.. code-block:: bash

  git clone https://gitlab.awi.de/tkleiner/grd2stream.git
  cd grd2stream
  ./bootsrap.sh
  # generic GMT *.grd file format
  ./configure --prefix=$HOME NETCDF_INC=/opt/local/include NETCDF_LIB=/opt/local/lib
  # all gdal readable file formats via GMT6 API (preferred)
  ./configure --prefix=$HOME --enable-gmt-api
  # check other compilers
  CC=clang ./configure --enable-gmt-api --enable-debug
  make && make dist


Contribute
----------

- Issue Tracker: https://github.com/tkleiner/grd2stream/issues
