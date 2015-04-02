#!/usr/bin/env python
##
## randowm surface
##

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import Rbf

# Generate data
np.random.seed(1981)
width, height = 300, 300
x, y, z = np.random.random((3,10))
x *= width
y *= height

#create a grid on which to interpolate data
xi, yi = np.mgrid[0:width:1j*width, 0:height:1j*height]

#interpolate the data with the matlab griddata function
interp = Rbf(x, y, z, function='linear')
zi = interp(xi, yi)


#create a matplotlib figure and adjust the width and heights
fig, ax = plt.subplots(subplot_kw=dict(frameon=True, xticks=[], yticks=[]))
ax.set_aspect('equal')

#create the contours and streamplot
CS = plt.contour(xi, yi, zi, linewidths=1, colors='b')
dy, dx = np.gradient(zi.T)
plt.streamplot(xi[:,0], yi[0,:], dx, dy, color='c', density=1, arrowsize=3)

plt.show()
