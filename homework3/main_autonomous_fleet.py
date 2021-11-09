import numpy as np
import matplotlib.pyplot as plt
from autonomous_fleet import AutonomousFleet

# Create a set of random points in two dimensions
# in the interval R = [-1, 1] x [-1, 1].
np.random.seed(42)
data = np.random.uniform(-1, 1, size=(100, 2))

# ------------------------------------------------------------------
# Create an AutonomousFleet object
obj = AutonomousFleet(data)

# ------------------------------------------------------------------
# 1. Nearest ships
# Find and report the s closest ships to the location using
# the Euclidean distance.
loc_point = data[2]
s = 4
nearest_ships = obj.nearest_ships(loc_point, 4)
print(
    f'{s} closest ships t the location of the ship {loc_point}: ')
for i in nearest_ships:
    print(i)

# ------------------------------------------------------------------
# 2. Avoiding collisions
# Find and return whether there are any ships within a square of
# arbitrary orientation of length r and center given by the
# location of a ship.
r = 0.5
ships_within_square = obj.ships_within_square(loc_point, r)
print(
    f'Ships within a square of length {r} and center \
given by the location of the ship {loc_point} rotated 90° \
respect to the center: ')
for i in ships_within_square:
    print(i)


# ------------------------------------------------------------------
# 3. Leading and lagging vessels
# Find and report all the ships that have locations
# with either maximal or minimal x- or y-coordinates.
leading_lagging_ships = obj.leading_lagging_ships()
print(
    'Ships that have locations with either maximal or minimal \
x- or y-coordinates: ')
print('N -> ', leading_lagging_ships[0])
print('S -> ', leading_lagging_ships[1])
print('E -> ', leading_lagging_ships[2])
print('W -> ', leading_lagging_ships[3])

# ------------------------------------------------------------------
# 4. Method that generates the k-d tree decomposition of the
# fleet with the location of each ship
obj.two_dimensional_map()
