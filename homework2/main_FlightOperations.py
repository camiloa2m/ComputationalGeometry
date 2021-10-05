import numpy as np
import pandas as pd
from FlightOperations import FlightOperations

texto = ' Corresponding test suite of class FlightOperations '
print(15*'#' + texto + 15*'#', end='\n\n')

# We read borders_CO.dat
borders_pts = np.fromfile('borders_CO.dat', sep=' ', dtype=float)
print('Data type borders_CO:', type(borders_pts))
print('Data shape borders_CO:', borders_pts.shape)

# Reshape borders_pts, column 0 : latitude, column 0 : longitude
borders_pts = np.reshape(borders_pts, (len(borders_pts)//2, 2))
print('New data shape borders_CO:', borders_pts.shape)

# We read airports_CO.dat
airports_data = pd.read_csv(
    'airports_CO.dat',
    sep="    ",
    names=['lat', 'lon', 'alt', 'city', 'dept', 'airport_name'],
    engine='python')
print('Data type airports_CO:', type(airports_data))
print('Data shape airports_CO:', airports_data.shape)
print('DataFrame airports_CO:\n', airports_data.head)

# We create an object FlightOper
FlightOper = FlightOperations(borders_pts, airports_data)

# Plot a figure of the airport and border locations from
# input data files and the 2D-voronoi diagram
FlightOper.plotVorMapAirports()
print()

# Find the airport names that cover the biggest and smallest areas
# of service and plot a figure of these areas.
airport_area_max, airport_area_min = FlightOper.minmaxCoveregeArea()
print('Airport that cover the biggest area of service:', airport_area_max)
print('Airport that cover the smallest area of service:', airport_area_min)
print()

# Find a point for a new airport and show a figure of
# the point and the coverage radius.
vertex_newAirport, rmax = FlightOper.newAiportArea(allcircles=False)
print('New airport location:', vertex_newAirport,
      '  Coverage radius:', rmax)
print()

# Find the Most and least crowded airports and show
# a figure of the location of this airports.
max_nb, min_nb = FlightOper.mostLeastCrowdedAirports()
print('Airport name:', max_nb[0],
      '  Number of neighboring airports:', max_nb[1])
print('Airport name:', min_nb[0],
      '  Number of neighboring airports:', min_nb[1])
print()

# Find airports that are too close to each other.
min_dist, mindist_points = FlightOper.closestPairOfAirports()
print('Closest pair airports position as a pair of locations: ',
      mindist_points[0], mindist_points[1])
print('Distance of the closest pair airports: ', min_dist)

print()
print(15*'#' + len(texto)*'#' + 15*'#')
