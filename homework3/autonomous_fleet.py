import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import KDTree


class AutonomousFleet:
    """Class that offers different options of autonomous ship
    sailing (motion) planning. The functionality evaluates
    aspects of autonomous ship navigation within a fleet.
    These aspects include global and local properties such as
    nearest and farthest ships, analysis of potential collisions,
    and drawing maps of the fleet.
    """

    def __init__(self, points):
        self.points = np.asarray(points)
        self.kd_tree = KDTree(points)

    def nearest_ships(self, loc_point, s):
        """Find the s closest ships to the loc_point using the
        Euclidean distance.

        Args:
            loc_point (list, tuple or array): point in the plane
            s (int): number of nearest neighbors to the location
            using the Euclidean distance.

        Returns:
            List: s closest ships to the loc_point using the Euclidean
            distance.
        """

        # if loc_point in points, s = s + 1
        dd, _ = self.kd_tree.query(loc_point, k=1)
        ss = s+1 if dd == 0.0 else s

        # dd: float or array of floats
        # The distances to the nearest neighbors
        # ii: integer or array of integers
        # The index of each neighbor
        dd, ii = self.kd_tree.query(loc_point, k=ss)

        plt.scatter(
            self.points[:, 0], self.points[:, 1],
            s=10, marker='D', c='gray')

        for i in ii:
            plt.scatter(
                self.points[i][0], self.points[i][1],
                s=10, marker='D')

        plt.scatter(
            loc_point[0], loc_point[1],
            s=25, marker='*', c='r', label=f'{loc_point}')

        plt.title(f'{s} closest ships to the location {loc_point}')
        plt.axis('equal')
        plt.legend()
        plt.show()

        return [self.points[i] for i in ii]

    def ships_within_square(self, loc_point, r):
        """Find whether there are any ships within a square of
        length r and center given by loc_point rotated pi/4
        respect to the center.

        Args:
            loc_point (list, tuple or array): point in the plane
            r (int): length of the square

        Returns:
            List:  list of the the neighbors of loc_point.
        """

        # Using the Minkowski norm with p = 1 we have a square
        # with sides parallel to the axes rotated pi/4.
        # List of the indices of the neighbors of loc_point
        querybp = self.kd_tree.query_ball_point(loc_point, r,  p=1)

        plt.scatter(
            self.points[:, 0], self.points[:, 1],
            s=10, marker='D', c='gray')

        for i in querybp:
            nearby_points = self.points[i]
            plt.scatter(
                nearby_points[0], nearby_points[1],
                s=10, marker='D')

        plt.scatter(
            loc_point[0], loc_point[1],
            s=25, marker='*', c='r', label=f'{loc_point}')

        plt.title(f'Ships within a square of length {r} and center\n\
given by the location of the ship {loc_point}.\nThe square is rotated $90°$ \
respect to the center')
        plt.axis('equal')
        plt.legend()
        plt.show()

        # if loc_point in points then loc_point in querybp_points
        querybp_points = [self.points[i] for i in querybp]

        return querybp_points

    def leading_lagging_ships(self):
        """Find and report all the ships that have locations with
        either maximal or minimal x- or y-coordinates.

        Returns:
            List: List of arrays. The List has length 4. [north points,
            south points, esat points, west points]
        """

        idxmaxY = np.where(self.points[:, 1] == self.kd_tree.maxes[1])
        idxminY = np.where(self.points[:, 1] == self.kd_tree.mins[1])
        idxmaxX = np.where(self.points[:, 0] == self.kd_tree.maxes[0])
        idxminX = np.where(self.points[:, 0] == self.kd_tree.mins[0])

        maxsY = idxmaxY[0]
        minsY = idxminY[0]
        maxsX = idxmaxX[0]
        minsX = idxminX[0]

        result = (
            self.points[maxsY, :],
            self.points[minsY, :],
            self.points[maxsX, :],
            self.points[minsX, :])

        plt.scatter(
            self.points[:, 0], self.points[:, 1],
            s=5, marker='D', c='gray')
        markss = '^v><'
        colorss = 'rrbb'
        for i in range(4):
            for j in range(len(result[i])):
                plt.scatter(
                    result[i][j][0],
                    result[i][j][1],
                    s=20, marker=markss[i], c=colorss[i])

        plt.title('Ships that have locations with either maximal\n\
or minimal x- or y-coordinates')
        plt.axis('equal')
        plt.show()

        return result

    class Node:
        """Node class"""

        def __init__(self, data):
            self.data = np.asarray(data)
            self.leftChild = None
            self.rightChhild = None
            self.split = ([], 0)  # (split node, split axis)
            # split axis = 0 : x, 1 : y

    def _kdtree(self, data, depth=0, mins=None, maxs=None, ax=None):
        """Private function that generate the kd tree unbalanced
        and make the figure of the-d tree decomposition in the
        object ax.
        """

        if data.size == 0:
            return None

        # Select axis based on depth
        axis = depth % 2

        if mins is None:
            offset = 0.25*np.array([1, 1])
            mins = data.min(0) - offset  # min element by axis x
        else:
            mins = np.asarray(mins)  # min element in axis = axis
        if maxs is None:
            offset = 0.25*np.array([1, 1])
            maxs = data.max(0) + offset  # max element by axis x
        else:
            maxs = np.asarray(maxs)  # max element in axis = axis

        # Sort point list and choose median as split point
        # axis = 0: sort by x, axis = 1: sort by y
        data = list(data)
        data.sort(key=lambda x: x[axis])

        data = np.asarray(data)

        median = len(data)//2  # choose median index

        # new maxs and mins
        maxs1 = maxs.copy()
        maxs1[axis] = data[median][axis]
        mins2 = mins.copy()
        mins2[axis] = data[median][axis]

        # plot splitting lines along the x(y)-coordinate.
        if depth % 2 == 0:
            c = 'orange'
            ax.plot([mins2[0], mins2[0]], [mins[1], maxs[1]], c=c)
        else:
            c = 'deepskyblue'
            ax.plot([mins[0], maxs[0]], [maxs1[1], maxs1[1]], c=c)

        # construct node
        node = self.Node(data)

        # (split node, axis is the the split axis)
        # axis values => 0 : x, 1 : y
        node.split = (data[median], axis)

        # construct children
        node.leftChild = self._kdtree(
            data[0:median], depth+1, mins, maxs1, ax=ax)
        node.rightChild = self._kdtree(
            data[median+1:], depth+1, mins2, maxs, ax=ax)

        return node

    def two_dimensional_map(self):
        """Generates a two-dimensional map with the location of
        each ship in the fleet.
        Make the k-d tree decompositionof the fleet. Dots
        represents ships, lines correspond to splitting lines
        along the x or y coordinate.
        """

        fig, ax = plt.subplots()

        ax.scatter(
            self.points[:, 0], self.points[:, 1],
            c='gray', s=10, marker='D',
            zorder=self.points.shape[0])

        tree = self._kdtree(self.points, ax=ax)

        ax.axis('equal')
        ax.plot()

        ax.set_title('k-d tree decomposition of the fleet\nDots \
represents ships, lines correspond to splitting lines\nalong the \
x or y coordinate')

        plt.show()


if __name__ == '__main__':
    print(30*'*', 'You are running the module', 30*'*')
    cwd = os.getcwd()  # Get the current working directory (cwd)
    print('Current working directory (cwd): ', cwd)
    files = os.listdir(cwd)  # Get all the files in that directory
    print("Files in cwd: %s" % (files))
