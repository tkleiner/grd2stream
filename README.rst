Implements a Runge-Kutta stream line integration based on gridded velocity
fields in native GMT grid format (vx.grd, vy.grd). The GMT5-API (or newer)
is not used at the moment.

.. code-block:: bash

    echo "0 0" | grd2stream vx.grd vy.grd | psxy -m -R -J ...


Example
-------

.. figure:: /docs/_static/map_FRIS_B13_line_source_small.png
   :scale: 50 %
   :alt: map of Filchner-Ronne Ice Shelf (FRIS) including ice core location B13

   Map of Filchner-Ronne Ice Shelf (FRIS) including ice core location B13. 
   Several stream lines as point or line sources from grd2stream and GMT mapping tools. 



ToDo:
-----
- Use GMT-API for file io!
- If GMT5 is installed with GDAL support, all GDAL file formats can be read (e.g., netCDF, GTiFF, etc.)


Installation
------------

Install grd2stream e.g. by

.. code-block:: bash

  ./configure --prefix=$HOME NETCDF_INC=/opt/local/include NETCDF_LIB=/opt/local/lib
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
  make && make dist
  


    

Contribute
----------

- Issue Tracker: https://gitlab.awi.de/tkleiner/grd2stream/issues
- Source Code: https://gitlab.awi.de/tkleiner/grd2stream/tree/master


