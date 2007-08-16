#!/bin/python
#
# integration example
#
# Written by Davide Bresolin on May, 2007
#

from ariadne import *
from math import *

vector=extract_vector
matrix=extract_matrix

# y' = z & z' = -y
dyn=AffineVectorField(matrix([[0,0,0],[0,0,1],[0,-1,0]]),vector([1,0,0]))

grid=Grid(Vector("[0.02,0.01,0.01]"))
block=LatticeBlock("[-2,352]x[-112,112]x[-112,112]")
fgrid=FiniteGrid(grid,block)
print fgrid

print "Creating reference set"
reference_set=GridMaskSet(fgrid)
for i in range(351):
	x = 0.02*i+0.001
	y = cos(x)
	pstr = "["+str(x)+","+str(x)+"]x["+str(y)+","+str(y)+"]x[0.001,0.001]"
	# print pstr
	point = Rectangle(pstr)
	reference_set.adjoin_outer_approximation(point)
	
print "reference_set.size(),capacity()=",reference_set.size(),reference_set.capacity()

print "Creating inital set"
initial_rectangle=Rectangle("[0.001,0.001]x[0.999,0.999]x[0.001,0.001]")
initial_list_set=RectangleListSet(3)
initial_grid_set=GridMaskSet(fgrid)
print "Adjoining initial rectangle"
initial_list_set.adjoin(initial_rectangle)
initial_grid_set.adjoin_outer_approximation(initial_rectangle)
print "initial_list_set.size()=",initial_list_set.size()
print "initial_grid_set.size(),capacity()=",initial_grid_set.size(),initial_grid_set.capacity()
chainreach_set=initial_grid_set

bounding_box=Rectangle("[-0.0,7.0]x[-1.1,1.1]x[-1.1,1.1]")

print "Creating bounding set"
bounding_set=GridMaskSet(fgrid)
print "Adjoining bounding rectangle"
bounding_set.adjoin_over_approximation(bounding_box)
print "bounding_set.size(),capacity()=",bounding_set.size(),bounding_set.capacity()

maximum_step_size=0.005;
lock_to_grid_time=10;
maximum_set_radius=0.1;

integrator=AffineIntegrator(maximum_step_size,lock_to_grid_time,maximum_set_radius)

print "Computing reach sets with Affine Integrator..."
reach_set=integrator.reach(dyn,initial_list_set,Rational(7))
print "reach_set.size()=",reach_set.size()
#hainreach_set=integrator.chainreach(dyn,initial_grid_set,bounding_set)
print "chainreach_set.size(),capacity()=",chainreach_set.size(),chainreach_set.capacity()

print "Exporting to postscript output...",
eps=EpsPlot()
eps.open("cos-affine.eps",bounding_box,0,1)

eps.set_line_style(True)
eps.set_fill_colour("white")
eps.write(bounding_box)
eps.set_fill_colour("green")
eps.write(chainreach_set)
eps.set_fill_colour("yellow")
eps.write(reach_set)
eps.set_fill_colour("magenta")
eps.write(reference_set)
eps.set_fill_colour("blue")
eps.write(initial_grid_set)
eps.close()

print "Done.\n"

#set_integrator_verbosity(6)
integrator=LohnerIntegrator(maximum_step_size,lock_to_grid_time,maximum_set_radius)

print "Computing reach sets with Lohner Integrator..."
reach_set=integrator.reach(dyn,initial_list_set,Rational(7))
print "reach_set.size()=",reach_set.size()
#chainreach_set=integrator.chainreach(dyn,initial_grid_set,bounding_set)
print "chainreach_set.size(),capacity()=",chainreach_set.size(),chainreach_set.capacity()

print "Exporting to postscript output...",
eps.open("cos-lohner.eps",bounding_box,0,1)

eps.set_line_style(True)
eps.set_fill_colour("white")
eps.write(bounding_box)
eps.set_fill_colour("green")
eps.write(chainreach_set)
eps.set_fill_colour("yellow")
eps.write(reach_set)
eps.set_fill_colour("magenta")
eps.write(reference_set)
eps.set_fill_colour("blue")
eps.write(initial_grid_set)
eps.close()

print "Done.\n"

#set_integrator_verbosity(6)
integrator=C1LohnerIntegrator(maximum_step_size,lock_to_grid_time,maximum_set_radius)

print "Computing reach sets with Lohner Integrator..."
reach_set=integrator.reach(dyn,initial_list_set,Rational(7))
print "reach_set.size()=",reach_set.size()
#chainreach_set=integrator.chainreach(dyn,initial_grid_set,bounding_set)
print "chainreach_set.size(),capacity()=",chainreach_set.size(),chainreach_set.capacity()

print "Exporting to postscript output...",
eps.open("cos-c1lohner.eps",bounding_box,0,1)

eps.set_line_style(True)
eps.set_fill_colour("white")
eps.write(bounding_box)
eps.set_fill_colour("green")
eps.write(chainreach_set)
eps.set_fill_colour("yellow")
eps.write(reach_set)
eps.set_fill_colour("magenta")
eps.write(reference_set)
eps.set_fill_colour("blue")
eps.write(initial_grid_set)
eps.close()

print "Done.\n"

set_integrator_verbosity(7)
maximum_step_size=0.05; lock_to_grid_time=1; maximum_set_radius=0.25

bounding_set=GridMaskSet(fgrid)
bounding_set.adjoin_over_approximation(Rectangle("[-0.0,2.0]x[-1.1,1.1]x[-1.1,1.1]"))

integrator=EulerIntegrator(maximum_step_size,lock_to_grid_time,maximum_set_radius)
print "Computing reach sets with Euler Integrator..."
reach_set=integrator.reach(dyn,initial_list_set,Rational(7))
print "reach_set.size()=",reach_set.size()
chainreach_set=integrator.chainreach(dyn,initial_grid_set,bounding_set)
print "chainreach_set.size(),capacity()=",chainreach_set.size(),chainreach_set.capacity()

print "Exporting to postscript output...",
eps.open("cos-euler.eps",bounding_box,0,1)

eps.set_line_style(True)
eps.set_fill_colour("white")
eps.write(bounding_box)
eps.set_fill_colour("green")
eps.write(reach_set)
eps.write(chainreach_set)
eps.set_fill_colour("magenta")
eps.write(reference_set)
eps.set_fill_colour("blue")
eps.write(initial_list_set)
eps.close()

print "Done.\n"



  	
