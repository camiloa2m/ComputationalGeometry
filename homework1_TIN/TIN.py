import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from numpy.core.arrayprint import printoptions
from scipy.spatial import Delaunay
from scipy.spatial import distance
from scipy import interpolate
from scipy.sparse import csgraph


class TIN:

    """
    n = number of points
    pts: numpyarray shape (n,2). It has n ordered pairs (x,y).
    elevs: numpyarray shape (n,1)
    """
    def __init__(self, pts, elevs):
        self.num_pts = pts.shape[0]
        self.pts = pts  # array of points (x,y). size = n x 2
        self.elevs = elevs
        self.delaunayTri = Delaunay(np.array(pts))

    def plotElevProfile3D(self):
        """Show a matplotlib's window with the respective plot
        """
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection='3d')
        ax.set_title("3D Elevation Profile")
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        X = self.pts[:, 0]
        Y = self.pts[:, 1]
        Z = self.elevs
        trisurf = ax.plot_trisurf(
            X, Y, Z,
            triangles=self.delaunayTri.simplices,
            cmap=cm.viridis)
        fig.colorbar(trisurf, shrink=0.5, aspect=3)
        plt.show()

    def linearInterpolation(self, x, y):
        """Return the linear interpolation of the point
        """
        # Indices of the points forming the simplices
        # in the triangulation. For 2-D,
        # the points are oriented counterclockwise.
        TIN_simplices = self.delaunayTri.simplices

        # Find the simplices containing the given points.
        # Indices of simplices containing each point.
        # Points outside the triangulation get the value -1.
        idx_simplice = self.delaunayTri.find_simplex((x, y))
        new_simplice = TIN_simplices[idx_simplice]

        X = (self.pts[new_simplice[0], 0],
             self.pts[new_simplice[1], 0],
             self.pts[new_simplice[2], 0])
        Y = (self.pts[new_simplice[0], 1],
             self.pts[new_simplice[1], 1],
             self.pts[new_simplice[2], 1])
        Z = (self.elevs[new_simplice[0]],
             self.elevs[new_simplice[1]],
             self.elevs[new_simplice[2]])

        f = interpolate.LinearNDInterpolator(list(zip(X, Y)), Z)
        z = f(x, y)

        return z

    def plotPointInProfile3D(self, x, y, z):
        """Show a matplotlib's window with the respective plot
        """
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection='3d')
        ax.set_title(f"Point in the 3D Elevation Profile")
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')

        X = self.pts[:, 0]
        Y = self.pts[:, 1]
        Z = self.elevs
        trisurf = ax.plot_trisurf(X, Y, Z,
                                  triangles=self.delaunayTri.simplices,
                                  cmap=cm.viridis)

        ax.scatter(x, y, z,
                   s=100, marker='d', c='r',
                   label=f'Point ({x}, {y}, {z})')
        ax.legend()
        fig.colorbar(trisurf, shrink=0.5, aspect=3)
        plt.show()

    def __areaTriangle(self, p0, p1, p2):
        """ Private method. Return the area of a tringle given its points
        """
        dst1 = distance.euclidean(p0, p1)
        dst2 = distance.euclidean(p1, p2)
        dst3 = distance.euclidean(p2, p0)
        s = (dst1 + dst2 + dst3) / 2  # calculate the semi-perimeter
        area = (s*(s-dst1)*(s-dst2)*(s-dst3)) ** 0.5  # calculate the area
        return area

    def largestAreaDrainageBasin(self, x, y):
        """Show a matplotlib's window with the largest area drainage
        basin given a planar location
        """
        # Indices of the points forming the simplices
        # in the triangulation. For 2-D,
        # the points are oriented counterclockwise.
        TIN_simplices = self.delaunayTri.simplices

        # Find the simplices containing the given points.
        # Indices of simplices containing each point.
        # Points outside the triangulation get the value -1.
        idx_simplice = self.delaunayTri.find_simplex((x, y))
        new_simplice = TIN_simplices[idx_simplice]

        # Indices of neighbor simplices for each simplex.
        # The kth neighbor is opposite to the kth vertex.
        # For simplices at the boundary, -1 denotes no neighbor.
        simplices_neighbors = self.delaunayTri.neighbors
        idx_simplice_neighbors = simplices_neighbors[idx_simplice]

        area_max = -1
        best_idx = 0
        for idx in idx_simplice_neighbors:
            if idx != -1:
                neighbor_simplice = TIN_simplices[idx]
                X = [self.pts[neighbor_simplice[0], 0],
                     self.pts[neighbor_simplice[1], 0],
                     self.pts[neighbor_simplice[2], 0]]
                Y = [self.pts[neighbor_simplice[0], 1],
                     self.pts[neighbor_simplice[1], 1],
                     self.pts[neighbor_simplice[2], 1]]
                Z = [self.elevs[neighbor_simplice[0]],
                     self.elevs[neighbor_simplice[1]],
                     self.elevs[neighbor_simplice[2]]]
                area = self.__areaTriangle(
                    self.pts[neighbor_simplice[0]],
                    self.pts[neighbor_simplice[1]],
                    self.elevs[neighbor_simplice[2]])
                if area > area_max:
                    area_max = area
                    best_idx = idx

        # First, plot the complite 3D elevation profile
        fig = plt.figure()
        ax = fig.add_subplot(1, 2, 1, projection='3d')
        lb1 = f"Largest area drainage basin given a planar \n\
location ({x}, {y}) in the 3D Elevation Profile "
        ax.set_title(lb1)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')

        X = self.pts[:, 0]
        Y = self.pts[:, 1]
        Z = self.elevs
        trisurf = ax.plot_trisurf(X, Y, Z,
                                  triangles=self.delaunayTri.simplices,
                                  cmap=cm.viridis)
        fig.colorbar(trisurf, shrink=0.5, aspect=3)

        # Second, plot the largest area drainage basin
        # given a planar location
        X = [self.pts[new_simplice[0], 0],
             self.pts[new_simplice[1], 0],
             self.pts[new_simplice[2], 0]]
        Y = [self.pts[new_simplice[0], 1],
             self.pts[new_simplice[1], 1],
             self.pts[new_simplice[2], 1]]
        Z = [self.elevs[new_simplice[0]],
             self.elevs[new_simplice[1]],
             self.elevs[new_simplice[2]]]
        ax.plot_trisurf(X, Y, Z, color='r')

        neighbor_simplice = TIN_simplices[best_idx]

        X = [self.pts[neighbor_simplice[0], 0],
             self.pts[neighbor_simplice[1], 0],
             self.pts[neighbor_simplice[2], 0]]
        Y = [self.pts[neighbor_simplice[0], 1],
             self.pts[neighbor_simplice[1], 1],
             self.pts[neighbor_simplice[2], 1]]
        Z = [self.elevs[neighbor_simplice[0]],
             self.elevs[neighbor_simplice[1]],
             self.elevs[neighbor_simplice[2]]]
        ax.plot_trisurf(X, Y, Z, color='r')

        # plot the largest area drainage basin given
        # a planar location in other subplot
        ax = fig.add_subplot(1, 2, 2, projection='3d')
        lb2 = f"Zoom in the Largest area\
drainage basin\n given a planar location ({x},{y})"
        ax.set_title(lb2)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_ylabel('z')

        X = [self.pts[new_simplice[0], 0],
             self.pts[new_simplice[1], 0],
             self.pts[new_simplice[2], 0]]
        Y = [self.pts[new_simplice[0], 1],
             self.pts[new_simplice[1], 1],
             self.pts[new_simplice[2], 1]]
        Z = [self.elevs[new_simplice[0]],
             self.elevs[new_simplice[1]],
             self.elevs[new_simplice[2]]]
        ax.plot_trisurf(X, Y, Z, color='r')

        neighbor_simplice = TIN_simplices[best_idx]

        X = [self.pts[neighbor_simplice[0], 0],
             self.pts[neighbor_simplice[1], 0],
             self.pts[neighbor_simplice[2], 0]]
        Y = [self.pts[neighbor_simplice[0], 1],
             self.pts[neighbor_simplice[1], 1],
             self.pts[neighbor_simplice[2], 1]]
        Z = [self.elevs[neighbor_simplice[0]],
             self.elevs[neighbor_simplice[1]],
             self.elevs[neighbor_simplice[2]]]
        ax.plot_trisurf(X, Y, Z, color='r')

        plt.show()

    def __calcAngle(self, a, b, c):
        """ Private method. Return the angle of the arc abc
        """
        ba = a - b
        bc = c - b
        cosine_angle = np.dot(ba, bc)/(np.linalg.norm(ba)*np.linalg.norm(bc))
        angle = np.arccos(cosine_angle)
        return angle

    def accuracyOfSampling(self):
        """ Return the accuracy of the sampling performed on the terrain
        """
        simplices = self.delaunayTri.simplices
        max_angle = 0
        min_angle = 180
        simp1 = []
        simp2 = []

        for simp in simplices:
            a, b, c = self.pts[simp[0]], self.pts[simp[1]], self.pts[simp[2]]
            angle_maxi = max([self.__calcAngle(a, b, c),
                              self.__calcAngle(b, c, a),
                              self.__calcAngle(c, a, b)])
            angle_mini = min([self.__calcAngle(a, b, c),
                              self.__calcAngle(b, c, a),
                              self.__calcAngle(c, a, b)])
            if angle_maxi > max_angle:
                max_angle = angle_maxi
                simp1 = simp
            if angle_mini < min_angle:
                min_angle = angle_mini
                simp2 = simp

        accu = min_angle/max_angle

        X = [self.pts[simp1[0], 0],
             self.pts[simp1[1], 0],
             self.pts[simp1[2], 0]]
        Y = [self.pts[simp1[0], 1],
             self.pts[simp1[1], 1],
             self.pts[simp1[2], 1]]
        plt.triplot(X, Y)

        X = [self.pts[simp2[0], 0],
             self.pts[simp2[1], 0],
             self.pts[simp2[2], 0]]
        Y = [self.pts[simp2[0], 1],
             self.pts[simp2[1], 1],
             self.pts[simp2[2], 1]]

        lb = 'Triangle with the largest and smallest\n\
angles in the triangulation'
        plt.title(lb)
        plt.triplot(X, Y, label=f'\nAccuracy of the sampling: {accu}')
        plt.legend(prop={'size': 8}, loc='upper right')
        plt.show()

        return accu

    def completePipeNetwork(self):
        """Show a matplotlib's window with the complete pipe network.
        Such network is defined as a planar graph connecting all the
        locations sampled; namely,
        it is the Euclidian minimum spanning tree of the sample points.
        """
        # We construct the dense, masked,and sparse representations
        # of the triangulation as follows
        graph_matrix = np.zeros((self.pts.shape[0], self.pts.shape[0]))
        for simp in self.delaunayTri.simplices:
            graph_matrix[simp[0], simp[1]] = distance.euclidean(
                self.pts[simp[0], :2],
                self.pts[simp[1], :2])
            graph_matrix[simp[1], simp[0]] = distance.euclidean(
                self.pts[simp[1], :2],
                self.pts[simp[0], :2])
            graph_matrix[simp[1], simp[2]] = distance.euclidean(
                self.pts[simp[1], :2],
                self.pts[simp[2], :2])
            graph_matrix[simp[2], simp[1]] = distance.euclidean(
                self.pts[simp[2], :2],
                self.pts[simp[1], :2])
            graph_matrix[simp[2], simp[0]] = distance.euclidean(
                self.pts[simp[2], :2],
                self.pts[simp[0], :2])
            graph_matrix[simp[0], simp[2]] = distance.euclidean(
                self.pts[simp[0], :2],
                self.pts[simp[2], :2])
        graph_masked = np.ma.masked_values(graph_matrix, 0)
        graph = csgraph.csgraph_from_dense(graph_matrix)

        # We make the minimum spanning tree of the undirected graph.
        # It's the N x N compressed-sparse representation (csr matrix) of the
        # undirected minimum spanning tree over the input
        mst = csgraph.minimum_spanning_tree(graph)

        # We make the adjacency matrix of mst
        adj_matrix_mst = mst.toarray().astype(np.float64)

        # We create the plot
        fig = plt.figure()
        # Draw a unstructured triangular grid as lines and/or markers.
        # The simplices attribute contains the indices of the points
        # in data that make up the triangle.
        plt.triplot(self.pts[:, 0],
                    self.pts[:, 1],
                    self.delaunayTri.simplices,
                    lw=0.5)
        plt.plot(self.pts[:, 0], self.pts[:, 1], 'mo', ms=2.5)
        for i in range(adj_matrix_mst.shape[0]):  # for rows
            for j in range(adj_matrix_mst.shape[0]):  # for colums
                if adj_matrix_mst[i][j] != 0:
                    plt.title(f"Complete pipe network")
                    plt.xlabel('x')
                    plt.ylabel('y')
                    X = self.pts[i, 0], self.pts[j, 0]
                    Y = self.pts[i, 1], self.pts[j, 1]
                    plt.plot(X, Y, 'm')
        plt.show()


if __name__ == '__main__':
    print(30*'*', 'You are running the module TIN!', 30*'*')
    cwd = os.getcwd()  # Get the current working directory (cwd)
    print('Current working directory (cwd): ', cwd)
    files = os.listdir(cwd)  # Get all the files in that directory
    print("Files in cwd: %s" % (files))
