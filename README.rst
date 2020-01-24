grd2stream

Implements a simple Runge-Kutta stream line integration based on gridded velocity
fields in native GMT grid format (vx.grd, vy.grd). The GMT5-API (or newer)
is not used at the moment.

.. code-block:: bash

    echo "0 0" | grd2stream vx.grd vy.grd | psxy -m -R -J ...


ToDo:
-----
- Use GMT-API for file io!
- If GMT5 is installed with GDAL support, all GDAL file formats can be read (e.g., netCDF, GTiFF, etc.)


Installation
------------

Install grd2stream e.g. by ::
  ./configure --prefix=$HOME ./configure NETCDF_INC=/opt/local/include NETCDF_LIB=/opt/local/lib
  make && make install

Packaging (Maintainer only)
---------------------------

Build grd2stream-X.X.X.tar.gz e.g. by ::
  git clone https://gitlab.awi.de/tkleiner/grd2stream.git
  cd grd2stream
  ./bootsrap.sh
  ./configure --prefix=$HOME ./configure NETCDF_INC=/opt/local/include NETCDF_LIB=/opt/local/lib
  make && make dist
  


    

Contribute
----------

- Issue Tracker: https://gitlab.awi.de/tkleiner/grd2stream/issues
- Source Code: https://gitlab.awi.de/tkleiner/grd2stream/master


