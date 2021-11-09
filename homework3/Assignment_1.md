## Introduction
So you've decided to spend your undergraduate internship at a robot motion planning company.
The company is widely known for their expertise in the control of the sailing of autonomous cargo
ships. These crewless vessels transport either containers or bulk cargo over navigable waters with
little or no human interaction. See Figure 1. Unfortunately, the last project on autonomous cargo
ships of the company has presented some problems related to the self-navigation system and help
is needed getting it back up again.
Given that you're a savvy intern with wide knowledge on algorithms and data structures, your
boss has assigned to you fixing the project problem(s). Solving the issues with the self-navigation
system of an autonomous cargo ship requires the implementation of a brand new interface. This
piece of software will be plugged into the rest of the software developed by other colleagues in the
company.

## Problem
Your job is then to provide a computational application that offers different options of autonomous
ship sailing (motion) planning. The functionality required evaluates aspects of autonomous ship
navigation within a fleet (a group of ships). These aspects include global and local properties such
as nearest and farthest ships, analysis of potential collisions, and drawing maps of the fleet.

## Test suite, with the following specifications:
We will treat locations of the ships on the map as defined by an ordered
pair \\[(x, y)\\] of coordinates and measure all distances using the standard Euclidean norm in the plane.
The methods must satisfy the following specifications:

1. Nearest ships. Communication between ships has a finite range. A ship has to know to which
ships send information relevant to the navigation route. Given the location of a vessel (point
in the plane), find and report the s closest ships to the location using the Euclidean distance.
If there is more than s closest to the location of the vessel, report only s of them, choosing
them randomly.
2. Avoiding collisions. As the ships navigate, they need to take into account if there are other
vessels within some distance so maneuvers can be performed and collisions can be avoided.
Find and return whether there are any ships within a square of arbitrary orientation of length
r and center given by the location of a ship. The orientation of the square must set by either
the library or the application. Only one range search needs to be performed.
3. Leading and lagging vessels. Suppose you want to explore areas in four directions: north (N),
south (S), west (W), and east (E). It would be advisable to send the ships that are either
lagging or leading within the 
eet to save fuel. Find and report all the ships that have locations
with either maximal or minimal x- or y-coordinates. Report all the ships if there is more than
one ship for each direction. It might happen that one ship can have maximal/minimal position
in more that one direction.
4. Implement a method that generates a two-dimensional map with the location of each ship in
the fleet. Represent ships as dots (or some other relevant symbol) and the splitting lines using
different colors for x- and y-coordinates.

The solution to each problem must generate a plot clearly showing and labeling the
answer; you may also print extra information to the terminal with further details. Regarding input
data, randomly generate around 10-100 points in the interval R = \\[[-1,1] \times [-1,1]\\].

## References
* J.J. Rincón, Universidad del Rosario. Midterm assignment 3. October 12, 2021. Computational and Differential Geometry 2021-2. Bogotá, Colombia
* M. de Berg, O. Cheong, M. van Kreveld, and M. Overmars. Computational Geometry: Algorithms
and Applications. Third edition. Springer-Verlag, 2008.
