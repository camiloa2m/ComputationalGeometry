import sys
from TIN import TIN
import numpy as np

# print(sys.path)

# We read the data
data = np.fromfile('pts1000c.dat', sep=' ', dtype=float)  # 1D-nparray
print('Data shape:', data.shape)

# we make a matrix (nx3) with data. n is the number of 3-tuples (x,y,z)
# col 0 = coord x, col 1 = coor y, col 2 = elevation
data = np.reshape(data, (len(data)//3, 3))
print('New data shape:', data.shape)

# We make tha array pts
pts = data[:, :2]
print('pts shape:', pts.shape)
# We make tha array elevs
elevs = data[:, 2]
print('elevs shape:', elevs.shape)

# We make a TIN object with the data
TIN_1 = TIN(pts, elevs)

# We make a 3D graph of the elevation profile
TIN_1.plotElevProfile3D()

# We make an interpolation function
# Using f , we make a elevation prediction
x, y = 9.319, 6.506
z = TIN_1.linearInterpolation(x, y)
print(f'Prediction of ({x}, {y}) = ', z)

# We plot the planar point whit its prediction in a
# 3D graph of the elevation profile
TIN_1.plotPointInProfile3D(x, y, z)

# Given a planar location (x,y) we find an plot its largest-area drainage basin
x = 0
y = 0
TIN_1.largestAreaDrainageBasin(x, y)

# We calculate the accuracy of the sampling performed on the terrain
accu = TIN_1.accuracyOfSampling()
print('Accuracy of the sampling = ', accu)

# We plot a complete pipe network.
# Such network is defined as a planar graph connecting all the
# locations sampled; # namely, it is the Euclidian minimum
# spanning tree of the sample points.
TIN_1.completePipeNetwork()
