## Introduction
The ministry of defense wants to build a server farm that is going to keep sensitive information. The server farm must fulfill a special requirement: there can be no traces of humidity or water
infiltration; otherwise data may be compromised. Therefore, the government hires a hydrology company to explore the construction site so they can locate all water and humidity sources; like
those coming from rain, for instance.
In order to map out the terrain the company has randomly sampled different parts of it, registering location and elevation. The final goal is to provide information on where to locate drain
pipes and design a drainage network so the farm can stay dry with no humidity issues that could damage the servers.
A drainage network keeps information on the hydrologic properties of the terrain; for example, water source locations and potential 
ood areas. The collection of water streams, sources, basins,
and sinks forms such a network, see Fig. 1.
![imagen](https://user-images.githubusercontent.com/33555617/136114409-f921ec8d-eaa0-452a-873e-bf0b77575e8e.png)
There is a problem, however. The hydrology company has no software to post-process the information and construct the whole drainage network. So they have decided to outsource this part
of the job to you! It is going to be your duty, then, to solve the problem of the drainage network.

## Problem
Design and implement an application, along with a test suite, that specifically addresses the functionality needed by the hydrology company. The functionality required evaluates
aspects of the hydrologic features of the construction site. Such features are local properties of drainage network such as points, lines, and regions of special interest to the 
ow of water due to gravity and obstacles present on the terrain

## Strategy
A good computational model is a digital elevation model (DEM), which is a computational model for representing terrain relief based on a finite number of points or samples. There are different
types of DEMs, but since the sample points on the terrain are not uniformly distributed the best DEM for the job is a triangulated irregular network.
A triangulated irregular network (TIN) is a polyhedral terrain; that is, the image of a piecewise-defined linear function with domain in two dimensions. A TIN is a triangulation of the location
points, as vertices of the triangulation, augmented with the elevation data. This data structure
gives rise to a piecewise-linear function that allows for interpolation of the elevations at any location (point), even though such location is not in the sample set.
Hence, a TIN is generated from the sampled points on the terrain by computing the Delaunay triangulation having vertices at the sample points. Thus, a TIN is constructed using algorithms for
Delaunay triangulations. Having this we can work on the implementation of the application.

## Test suite, with the following specifications:
1. Plot the elevation profile using different colors for different ranges of elevations. Given the set of sample points on the terrain, plot its elevation profile as a three dimensional graph.
1. Given the planar location of a point and its corresponding elevation value. Do this by means of a linear interpolation using the sample data.
1. Given the planar location of a point and its largest-area drainage basin. The drainage basin is a quadrilateral, composed of two sampled triangles, that contains the planar location.
1. Measure the accuracy of the sampling performed on the terrain. That is, find the largest and smallest angles in the triangulation and report their values and ratio.
1. Find a complete pipe network. Such network is defined as a planar graph connecting all the locations sampled; namely, it is the Euclidian minimum spanning tree of the sample points.

## References
* J.J. Rincón, Universidad del Rosario. Midterm assignment 1. September 3, 2021. Computational and Differential Geometry 2021-2. Bogotá, Colombia.
* L. De Floriani, P. Magillo, and E. Puppo. Applications of Computational Geometry to Geographic Information Systems. In J.-R. Sack and J. Urrutia, editors, Handbook of Computational Geometry,
chapter 7, pages 333-388. North Holland, 2000.
* L. Floriani and P. Magillo. Algorithms for Visibility Computation on Terrains: A Survey. Environment and Planning B: Planning and Design 30, issue 5, pages 709-728 (2003).
