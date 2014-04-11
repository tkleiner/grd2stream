#!/usr/bin/env python
##
## real data
##

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import Rbf
from optparse import OptionParser  # command line options

try:
    from netCDF4 import Dataset as CDF
except:
    from netCDF3 import Dataset as CDF




## Set up the option parser
desc = """Bla Bla...."""

parser = OptionParser(description=desc)
parser.usage = "usage: %prog [options] in.nc out.nc"
# parser.add_option("-f","--fill_nan", 
#                   help="fill missing data (e.g. ocean grid cells)",
#                   dest="fill_nan",
#                   default=False, 
#                   action="store_true")

# parser.add_option("-K","--Kelvin", 
#                   help="output in kelvin units",
#                   dest="kelvin",
#                   default=False, 
#                   action="store_true")

(opts, args) = parser.parse_args()

if len(args) == 1:
    ncfile = args[0]
else:
    print('wrong number arguments, 1 expected')
    parser.print_help()
    exit(0)


if __name__ == "__main__": 

    ####################################################################
    ## get data from ncfile
    ####################################################################

    ## open netCDF file in append mode
    nc = CDF(ncfile,'r')


    ## check for data
    names = ['u', 'velx', 'vel_x', 'vx']
    try:
        ## go through variables and look for 'grid_mapping' attribute
        for varname in names:
            if varname in nc.variables.keys():
                print('Found information in variable "%s", using it' % varname)
                u_var = nc.variables[varname]
                u = u_var[:]
                uname = varname
                exit
    except:
        print('Error: No information found')
        exit(-1)


    ## check for data
    names = ['v', 'vely', 'vel_y', 'vy']
    try:
        ## go through variables and look for 'grid_mapping' attribute
        for varname in names:
            if varname in nc.variables.keys():
                print('Found information in variable "%s", using it' % varname)
                v_var = nc.variables[varname]
                v = v_var[:]
                vname = varname
                exit
    except:
        print('Error: No information found')
        exit(-1)



    ## read dimensions
    ndims = len(u_var.shape)
    if ( ndims == 4):
        print ('Assuming time as first dimension')
        
        tdim = u_var.dimensions[0]
        zdim = u_var.dimensions[1]
        ydim = u_var.dimensions[2]
        xdim = u_var.dimensions[3]
        
        
        ## coordinate variable in time-direction
        time_var = nc.variables[tdim]
        time = time_var[:]
        nt = len(time)
        
        ## coordinate variable in x-direction
        x_var = nc.variables[xdim]
        x = x_var[:]
        nx = len(x) 
        
        ## coordinate variable in y-direction
        y_var = nc.variables[ydim]
        y = y_var[:]
        ny = len(y)
    
        ## coordinate variable in y-direction
        z_var = nc.variables[zdim]
        z = z_var[:]
        nz = len(z)
        
        
    else:
        print ('Error')
        exit(-1)
        
    vx = np.squeeze(u[nt-1,nz-1,:,:])*(3600*24*365.25)
    vy = np.squeeze(v[nt-1,nz-1,:,:])*(3600*24*365.25)


    nc.close()

    mag=np.sqrt(vx*vx+vy*vy)
    vx_x,vx_y = np.gradient(vx)
    vy_x,vy_y = np.gradient(vy)
    div=np.sqrt(vx_x*vx_x+vy_y*vy_y)
    
    #create a matplotlib figure and adjust the width and heights
    fig, ax = plt.subplots(subplot_kw=dict(frameon=True, xticks=[], yticks=[]))
    ax.set_aspect('equal')

    #create the contours and streamplot
    CS = plt.contour(x, y, mag, linewidths=1, colors='b')
    DS = plt.contour(x, y, div, linewidths=1, colors='r')
    plt.streamplot(x, y, vx, vy, color='c', density=3, arrowsize=1)
    plt.show()

