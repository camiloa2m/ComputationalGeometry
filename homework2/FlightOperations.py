import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d
from scipy.spatial import distance
from scipy.spatial import ConvexHull
from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection


class FlightOperations:

    """
    borders_pts: points of the country borders (lat, lon),  size (numpoints,2)
    airports_data: aiports data, pandasDataframe with columns
    ['lat', 'lon', 'alt', 'city', 'dept', 'airport_name']
    """
    def __init__(self, borders_pts, airports_data):
        # num_airports: number of airports
        self.num_airports = airports_data.shape[0]
        # aiports data, pandasDataframe
        self.airports_data = airports_data
        # points of the country borders (lat, lon)
        self.borders_pts = borders_pts
        # points of the airports (lon, lat)
        self.pts_airports = self._make_pts_airports(airports_data)
        # Voronoi diagram of pts_airports
        self.vor = Voronoi(self.pts_airports)
        # convex hull
        self.hull = ConvexHull(self.pts_airports)

    def _make_pts_airports(self, airports_data):
        '''Private method
        Return an array of the airports's points'''
        pts_airports = []
        for i in range(airports_data.shape[0]):
            pts_airports.append([airports_data.iloc[i, 1],
                                 airports_data.iloc[i, 0]])
        pts_airports = np.array(pts_airports)
        return pts_airports

    def plotVorMapAirports(self):
        '''Show a figure of the airport and border locations from
        input data files and the 2D-voronoi diagram'''
        # Plot the 2D-voronoi diagram
        fig = voronoi_plot_2d(
            self.vor, show_vertices=False,
            point_size=1,
            line_colors=['orange', 'limegreen', 'royalblue',
                         'magenta', 'darkviolet', 'teal'],
            zorder=2)

        # Plot the border
        plt.plot(self.borders_pts[:, 1],
                 self.borders_pts[:, 0],
                 '#808080', zorder=1)

        # Plot aiport points in a diferent color
        plt.scatter(self.airports_data['lon'],
                    self.airports_data['lat'],
                    s=3, c='k', zorder=3)

        title0 = 'Colombia\'s two-dimensional Voronoi diagram\n\
for flight space with airports as sites'

        plt.title(title0)
        plt.xlabel('Longitude (degrees)')
        plt.ylabel('Latitude (degrees)')
        plt.axis('equal')
        plt.grid()
        plt.xlim(-85, -65)
        plt.ylim(-6, 18)
        plt.show()

    def _polyArea(self, x, y):
        '''Private method
        Determine the area of a simple polygon using Shoelace formula'''
        area = 0.5*np.abs(np.dot(x, np.roll(y, 1))-np.dot(y, np.roll(x, 1)))
        return area

    def minmaxCoveregeArea(self):
        '''Find the airport names that cover the biggest and smallest
        areas of service and show a figure of these areas.
        Return:
        (airport_area_max_name, airport_area_min_name)'''
        max_area = 0
        max_idreg = 0
        min_area = 1e10
        min_idreg = 0
        for idreg in self.vor.point_region:
            region = self.vor.regions[idreg]
            if -1 not in region:
                # polygon defined by points
                x = [self.vor.vertices[k][0] for k in region]
                y = [self.vor.vertices[k][1] for k in region]
                area = self._polyArea(x, y)
                if area > max_area:
                    max_area = area
                    max_idreg = idreg
                if area < min_area:
                    min_area = area
                    min_idreg = idreg

        # Plot the 2D-voronoi diagram
        fig = voronoi_plot_2d(
            self.vor, show_vertices=False,
            point_size=1,
            line_colors=['orange', 'limegreen', 'royalblue',
                         'magenta', 'darkviolet', 'teal'],
            zorder=2)

        # Plot the border
        plt.plot(self.borders_pts[:, 1],
                 self.borders_pts[:, 0],
                 '#808080', zorder=1)

        # Plot aiport points in a diferent color
        plt.scatter(self.airports_data['lon'],
                    self.airports_data['lat'],
                    s=3, c='k', zorder=3)

        idxmax = np.where(self.vor.point_region == max_idreg)[0][0]
        idxmin = np.where(self.vor.point_region == min_idreg)[0][0]
        airport_area_max = self.airports_data.iloc[idxmax]
        airport_area_min = self.airports_data.iloc[idxmin]

        # Plot the biggest aiport area of service
        x = [self.vor.vertices[k][0] for k in self.vor.regions[max_idreg]]
        y = [self.vor.vertices[k][1] for k in self.vor.regions[max_idreg]]
        lb1 = f'Biggest aiport area of service: {airport_area_max[5]}'
        plt.fill(x, y, label=lb1)
        plt.legend()

        # Plot the smallest aiport area of service
        x = [self.vor.vertices[k][0] for k in self.vor.regions[min_idreg]]
        y = [self.vor.vertices[k][1] for k in self.vor.regions[min_idreg]]
        lb2 = f'Smallest aiport area of service: {airport_area_min[5]}'
        plt.fill(x, y, label=lb2)
        plt.legend()

        plt.title('Airports with the biggest and smallest areas of service')
        plt.xlabel('Longitude (degrees)')
        plt.ylabel('Latitude (degrees)')
        plt.axis('equal')
        plt.grid()
        plt.xlim(-85, -65)
        plt.ylim(-6, 18)
        plt.show()

        return airport_area_max[5], airport_area_min[5]  # airport names

    def _makeLeftTurn(self, p0, p1, p2):
        '''Private method
        Return True if a point p1 is left of p2 with respect to p0, that is
        vector (p0,p2) turn left of vector (p0,p1), otherwise False'''
        det = (p1[0]-p0[0])*(p2[1]-p0[1])-(p2[0]-p0[0])*(p1[1]-p0[1])
        izq = det >= 0
        return izq

    def _insideCH(self, pt):
        '''Private method
        Return True if a point pt is inside or on the convex hull
        otherwise False'''
        hull_vert = self.hull.vertices
        hull_vertRoll = np.roll(self.hull.vertices, -1)
        p2 = pt
        inside_CH = True
        for i in range(len(hull_vert)):
            p0 = self.vor.points[hull_vert[i]]
            p1 = self.pts_airports[hull_vertRoll[i]]
            if not self._makeLeftTurn(p0, p1, p2):
                inside_CH = False
                return inside_CH
        return inside_CH

    def newAiportArea(self, allcircles=False):
        ''' Find a point for a new airport and show
        a figure of the point and the coverage radius.
        Given a set of points (airports) in the plane, find a largest circle
        centered within their convex hull and enclosing none of them.
        Return:
        (Point_of_the_new_airport, coverage_radius)'''
        patches = []
        rmax = -1
        idp_rmax = -1
        idv_rmax = -1
        for i in range(1, len(self.vor.ridge_points)-1):
            idp = self.vor.ridge_points[i][1]
            for idv in self.vor.ridge_vertices[i]:  # iters 2 times
                if idv != -1 and self._insideCH(self.vor.vertices[idv]):
                    radius = distance.euclidean(
                        self.vor.vertices[idv], self.vor.points[idp])
                    if radius > rmax:
                        rmax = radius
                        idp_rmax = idp
                        idv_rmax = idv
                    if allcircles:
                        circle = Circle(self.vor.vertices[idv],
                                        radius=radius,
                                        ec="none")
                        patches.append(circle)

        new_vertex = self.vor.vertices[idv_rmax]

        fig, ax = plt.subplots()

        if allcircles:
            colors = np.linspace(0, 1, len(patches))
            collection = PatchCollection(
                patches, cmap=plt.cm.hsv, alpha=0.1, zorder=2)
            collection.set_array(colors)
            # Add circles collection
            ax.add_collection(collection)

        # Add circle of max radius
        circle = Circle(
            new_vertex, radius=rmax, ec="none", alpha=0.3, zorder=4)
        ax.add_patch(circle)

        # Plot the border
        ax.plot(self.borders_pts[:, 1],
                self.borders_pts[0:, 0],
                '#808080', zorder=1)

        # Plot aiport points in a diferent color
        ax.scatter(self.airports_data['lon'],
                   self.airports_data['lat'],
                   s=3, c='k', zorder=3)

        # Plot new aiport point
        ax.scatter(new_vertex[0],
                   new_vertex[1],
                   s=10, c='r', zorder=4,
                   label=f'Point: {new_vertex}\nRadius: {rmax}')
        ax.legend()

        ax.set_xlabel('Longitude (degrees)')
        ax.set_ylabel('Latitude (degrees)')
        plt.title('New airport and its coverage radius')
        plt.axis('equal')
        plt.grid()
        plt.xlim(-85, -65)
        plt.ylim(-6, 18)
        plt.show()

        return self.vor.vertices[idv_rmax], rmax

    def mostLeastCrowdedAirports(self):
        '''Find the Most and least crowded airports and show
        a figure of the location of this airports.
        Return:
        ((Most_crowded_airport_name, Number_of_neighboring_airports),
        (Least_crowded_airport_name, Number_of_neighboring_airports))'''
        neighbors = {}  # key:idxpoint, value:list of idx neighbors
        for k in range(len(self.pts_airports)):
            neighbors[k] = []

        for r_pts in self.vor.ridge_points:
            pt1, pt2 = r_pts
            neighbors[pt1].append(pt2)
            neighbors[pt2].append(pt1)

        max_nb = -1
        id_pt_maxnb = -1
        min_nb = 1e6
        id_pt_minnb = -1
        for k in neighbors.keys():
            num_nb = len(neighbors[k])
            if num_nb > max_nb:
                max_nb = num_nb
                id_pt_maxnb = k
            if num_nb < min_nb:
                min_nb = num_nb
                id_pt_minnb = k

        # Plot the 2D-voronoi diagram
        fig = voronoi_plot_2d(
            self.vor, show_vertices=False,
            point_size=1,
            line_colors=['orange', 'limegreen', 'royalblue',
                         'magenta', 'darkviolet', 'teal'],
            zorder=2)

        # Plot the border
        plt.plot(self.borders_pts[:, 1],
                 self.borders_pts[:, 0],
                 '#808080', zorder=1)

        # Plot aiport points in a diferent color
        plt.scatter(self.airports_data['lon'],
                    self.airports_data['lat'],
                    s=3, c='k', zorder=3)

        airport_maxnb = self.airports_data.iloc[id_pt_maxnb]  # aiport data
        airport_minnb = self.airports_data.iloc[id_pt_minnb]  # aiport data

        # Plot most crowded airport point.
        x = (self.vor.points[id_pt_maxnb][0])
        y = (self.vor.points[id_pt_maxnb][1])
        lb1 = f'Most crowded airport: {airport_maxnb[5]}'
        lb1 = lb1 + f'\nNumber of neighboring airports: {max_nb}'
        plt.scatter(x, y, s=10, c='r', zorder=4, label=lb1)
        plt.legend()

        # Plot least crowded airport point.
        x = (self.vor.points[id_pt_minnb][0])
        y = (self.vor.points[id_pt_minnb][1])
        lb2 = f'Least crowded airport: {airport_minnb[5]}'
        lb2 = lb2 + f'\nNumber of neighboring airports: {min_nb}'
        plt.scatter(x, y, s=10, c='c', zorder=4, label=lb2)
        plt.legend()

        plt.title('Most and least crowded airports')
        plt.xlabel('Longitude (degrees)')
        plt.ylabel('Latitude (degrees)')
        plt.axis('equal')
        plt.grid()
        plt.xlim(-85, -65)
        plt.ylim(-6, 18)
        plt.show()

        return (airport_maxnb[5], max_nb), (airport_minnb[5], min_nb)

    def closestPairOfAirports(self):
        '''Find airports that are too close to each other.
        Return: (min_dist, [pt1, pt2])
        pt1 and pt2 are the two close airports (in coordinates)
        '''
        # Closest pair of points problem
        min_dist = 1e10
        mindist_idpoints = []
        idpt1 = -1
        idpt2 = -1
        for r_pts in self.vor.ridge_points:
            idpt1, idpt2 = r_pts
            pt1, pt2 = self.pts_airports[idpt1], self.pts_airports[idpt2]
            dist = distance.euclidean(pt1, pt2)
            if dist < min_dist:
                min_dist = dist
                mindist_idpoints = [idpt1, idpt2]

        # Plot the 2D-voronoi diagram
        fig = voronoi_plot_2d(
            self.vor, show_vertices=False,
            point_size=1,
            line_colors=['orange', 'limegreen', 'royalblue',
                         'magenta', 'darkviolet', 'teal'],
            zorder=2)

        # Plot the border
        plt.plot(self.borders_pts[:, 1],
                 self.borders_pts[:, 0],
                 '#808080', zorder=1)

        # Plot aiport points in a diferent color
        plt.scatter(self.airports_data['lon'],
                    self.airports_data['lat'],
                    s=3, c='k', zorder=3)

        pt1 = self.pts_airports[mindist_idpoints[0]]
        pt2 = self.pts_airports[mindist_idpoints[1]]
        pt1_data = self.airports_data.iloc[mindist_idpoints[0]]  # aiport data
        pt2_data = self.airports_data.iloc[mindist_idpoints[1]]  # aiport data

        # Plot closest pair of airports
        plt.scatter(pt1[0], pt1[1],
                    s=10, c='r', zorder=4, label=f'{pt1_data[5]}: {pt1}')
        plt.scatter(pt2[0], pt2[1],
                    s=10, c='c', zorder=4, label=f'{pt2_data[5]}: {pt2}')
        plt.legend()

        plt.xlabel('Longitude (degrees)')
        plt.ylabel('Latitude (degrees)')
        plt.title('Closest pair of airports')
        plt.axis('equal')
        plt.grid()
        plt.xlim(-85, -65)
        plt.ylim(-6, 18)
        plt.show()

        return min_dist, [pt1, pt2]


if __name__ == '__main__':
    print(30*'*', 'You are running the module', 30*'*')
    cwd = os.getcwd()  # Get the current working directory (cwd)
    print('Current working directory (cwd): ', cwd)
    files = os.listdir(cwd)  # Get all the files in that directory
    print("Files in cwd: %s" % (files))
