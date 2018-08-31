#!/usr/bin/env python3
"""
Purpose: Example of locating earthquakes in a 2D medium with xcloc.
Copyright: Ben Baker distributed under the MIT license.
"""
from numpy.random import rand
from numpy.random import seed
from numpy import sqrt
from numpy import meshgrid
from numpy import asarray
from numpy import copy
from numpy import linspace
import sys
from acoustic2D import Acoustic2D
#sys.path.insert(0, './..') # This can be used if your PYTHONPATH isn't set
from xcloc import xcloc
from xclocTypes import xclocTypes as xctypes
import matplotlib.pyplot as plt


# Model parameters
nrec = 20       # Number of receivers
nsrc = 2        # Number of sources
dt = 1.0/6000.0 # Sampling rate is 6000 Hz
fcent = 400.0   # Dominant resolution is vmin/fcent
x0 = 0.0        # Model origin in x is 0 km
x1 = 1000.0     # Model extent in x is 1 km
y0 = 0.0        # Model origin in y is 0 km
y1 = 1000.0     # Model origin in z is 1 km
z0 = 0.0        # This is 2d
z1 = 0.0        # This is 2d
nx = 512        # Number of grid points in x
ny = 512        # Number of grid points in y
nz = 1          # Number of grid points in z
ngrd = nx*ny*nz # Number of grid points in model
vel = 3100.0    # P velocity (m/s)
rho = 2700.0    # Density (kg/m**3)
Q = 9.e2        # Qp
tmodel = 0.5    # Max modeling time is the travel time from the furthest
                # point in the medium to the receiver plus some.
srcScale = asarray([1., 1.1]) # Magnitude of each source
npts = int(round(tmodel/dt)) + 1 # Number of points in signals
dx = (x1 - x0)/(nx - 1) 
dy = (y1 - y0)/(ny - 1)
dz = 0.0
# Randomly space the receivers in the 2d model
seed(4343)
xrec = rand(nrec, 3) # Receiver positions in [0,1)
xrec[:,0] = x0 + (x1 - x0)*xrec[:,0]
xrec[:,1] = y0 + (y1 - y0)*xrec[:,1]
xrec[:,2] = z0 + (z1 - z0)*xrec[:,2]
# Randomly place source in the model
xs = rand(nsrc, 3)   # Source positions in [0,1)
xs[:,0] = x0 + (x1 - x0)*xs[:,0]
xs[:,1] = y0 + (y1 - y0)*xs[:,1]
xs[:,2] = z0 + (z1 - z0)*xs[:,2]
# Tabulate the source indices
srcIndex = []
for isrc in range(nsrc):
    isx = int((xs[isrc,0] - x0)/dx + 0.5)
    isy = int((xs[isrc,1] - y0)/dy + 0.5)
    isz = 0
    if (abs(dz) > 0):
        isz = int((xs[isrc,2] - z0)/dz + 0.5)
        sys.exit("Plotting will fail in full 3d") #print("This isn't done right")
    # Will need [nx x ny] matrix
    srcIndex.append(isx*ny + isy)
#for i in range(nrec):
print(xs)
print(xrec)
# Initialize xcloc
print("Initializing xcloc...")
nsignals = nrec  
xcloc = xcloc(xcloc_path=['/home/bakerb25/C/xcloc/lib'])
xcloc.initialize(dt, npts, ngrd, nsignals=nsignals)
# Create 2D Acoustic Green's functions assuming a line source in z
ac2d = Acoustic2D(vel=vel, rho=rho, Q=Q)
ricker = ac2d.ricker(npts, dt, fcent)
Gmat = None
for i in range(nsrc):
    print("Computing Green's functions for source: %d"%(i+1))
    G = ac2d.computeGreensLineSource(npts, dt, xs[i,:], xrec, stf=ricker)
    G = G*srcScale[i]
    if (Gmat is None):
        Gmat = copy(G)
    else:
        Gmat = Gmat + G
print("Setting signals...")
xcloc.setSignals(Gmat)
# Compute the reciprocity based travel time fields
print("Setting travel time tables...")
x = linspace(x0, x1, nx)
y = linspace(y0, y1, ny)
xv, yv = meshgrid(x, y, indexing='ij') # [nx, ny] matrix
for i in range(nrec):
    dmat = sqrt( (xrec[i,0] - xv)**2 + (xrec[i,1] - yv)**2 )
    t = dmat.flatten(order='C')/vel
    xcloc.setTable(i, t) 
print("Computing DSM image...")
xcloc.compute()
print("Getting DSM image...")
maxIndex, maxValue = xcloc.getImageMax()
if (not maxIndex in srcIndex):
    sys.exit("Failed to locate source") #print("louckyy")
print(maxIndex, maxValue, srcIndex)
image = xcloc.getImage()
image = image.reshape(xv.shape, order='C')
# Destroy xcloc
print("Destroying xcloc...")
xcloc.finalize()
# Plot it
image = image.reshape(xv.shape, order='C')
plt.imshow(image.T, cmap=plt.viridis(),
           extent=(min(x), max(x), max(y), min(y)))
plt.xticks(linspace(min(x), max(x), 11))
plt.yticks(linspace(max(y), min(y), 11))
plt.xlim([min(x), max(x)])
plt.ylim([max(y), min(x)])
#plt.colorbar()
for i in range(nrec):
    plt.scatter(xrec[i,0], xrec[i,1], marker='v', c='red', s=100)
plt.xlabel('Offset')
plt.ylabel("Depth")
plt.title("DSM Image")
plt.show()
