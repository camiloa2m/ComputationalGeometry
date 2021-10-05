## Introduction
An essential part in aviation is the planning of flights. Given a departure and a destination, a flight plan is a carefully designed route to take safely an aircraft between those two points, within an estimated time en route. As expected, not every flight might be pleasant so a flight plan must include contingency procedures, in case an emergency ensues. For that, several pieces of information
are included: estimated time en route, alternate airports in case of detouring, type of ight, crew's information, number of passengers on board, and aircraft specifications. Flight planning is extremely important, especially when ying over unexplored areas, as they provide vital information in case of mechanical failure, bad weather, or when rescuing operations are needed.
Also known as ight path, a route consists of line segments located on a map that aid the aircraft to navigate from the departure to the destination airports. Although there are several types of routes, with diferent attributes such as altitude and airspeed, we will consider a route to be defined as a sequence of two-dimensional line segments on a map. A single flight plan may include more than one ight path; that is, diferent routes for diferent scenarios.
Flight planning also involves the identification of airports which can be own to in case of unexpected conditions at the destination airport or along the route. The flight plan should only include alternate airports or locations which can be reached with the current fuel load and that have the capabilities necessary to handle the type of aircraft being own.

## Problem
Design and implement a computational application that addresses features of a flight plan. The functionality required evaluates aspects of the 
ight plan, given a set of airport sites, and departure and destination points of the flight. These aspects include global and local properties
such as preprocessing of airport sites, nearest and farthest airports, types of flight paths or routes, potential threats, and desolated areas.
We will treat locations or points on the map as defined by an ordered pair (latitude and longitude) and measure all distances using the standard
Euclidean norm in the plane. The methods must satisfy the following specifications:
1. Airport coverage area. It is always important to know how large and what population is
serviced by a given airport; this helps minimizing the ration of fuel to passenger number.
Find the airport names that cover the biggest and smallest areas of service. Notice that this
can only be computed for airports with finite area of service.
1. Build new airports. Suppose you want to explore desolated areas to gather information on
where it will be beneficial to build new airports, regardless of population and other factors.
You can find such information by finding the largest circle centered within the convex hull
and enclosing none of the already existing airports. Report this circle as the location of its
center and corresponding radius.
1. Most and least crowded airports. In case of an emergency landing it is important to know
how populated are the airports in order to cause the less damage. Report the names, not
locations, of the airports with the smallest and the largest number of neighboring airports.
Make sure that this search problem is solved in linear time in the structural complexity of the
Voronoi diagram.
