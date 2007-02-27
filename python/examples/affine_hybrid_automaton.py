#!/usr/bin/python

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

from ariadne import *

print dir()

space=PolyhedralSet(Rectangle("[-7.5,7.5]x[-7.5,7.5]"))
identity=AffineMap(Matrix("[1,0;0,1]"),Vector("[0,0]"))
dynamic=AffineVectorField(Matrix("[-2,-1; 1,-2]"),Vector("[-1,0]"))
jump_set = PolyhedralSet(Rectangle("[-7.5,7.5]x[-3,-2]"))
reset = AffineMap(Matrix("[1,0;0,-1]"),Vector("[0,0]"))

automaton=HybridAutomaton("Affine hybrid automaton")
mode_id=2
mode=automaton.new_mode(mode_id,dynamic,space)
event_id=3
transition=automaton.new_transition(event_id,mode_id,mode_id,reset,jump_set)
print automaton
#print automaton.invariant()

initial_rectangle=Rectangle("[-6.96875,-6.9375]x[-6.96875,-6.9375]");
bounding_box=Rectangle("[-8,8]x[-8,8]")


grid=RegularGrid(Vector("[0.125,0.125]"))
block=LatticeBlock("[-64,64]x[-64,64]")
fgrid=FiniteGrid(grid,block)
print fgrid

print "Creating initial hybrid set",
initial_set=HybridGridMaskSet()
initial_set.new_location(mode_id,fgrid)
initial_set[mode_id].adjoin_over_approximation(initial_rectangle)
print "initial_set.locations() =",initial_set.locations()
print initial_set[mode_id].size(),initial_set[mode_id].capacity()

print "Creating initial hybrid cell list set",
initial_cell_list_set=HybridGridCellListSet()
initial_cell_list_set.new_location(mode_id,grid)
initial_cell_list_set[mode_id].adjoin_over_approximation(initial_rectangle)
print initial_cell_list_set[mode_id].size()

print "Creating bounding hybrid set",
bounding_set=HybridGridMaskSet();
bounding_set.new_location(mode_id,fgrid)
bounding_set[mode_id].adjoin_over_approximation(bounding_box);
print bounding_set[mode_id].size(),bounding_set[mode_id].capacity()

apply=Applicator()
#integrator=LohnerIntegrator(0.125,0.5,0.25);
integrator=AffineIntegrator(0.125,0.5,0.25);
hybrid_evolver=HybridEvolver(apply,integrator);

print "Computing single discrete step"
discrete_step_set=hybrid_evolver.discrete_step(automaton,initial_cell_list_set)
print "Computing continuous chainreach set"
continuous_chainreach_set=hybrid_evolver.continuous_chainreach(automaton,initial_set,bounding_set)
print "Computing chainreach set"
chainreach_set=hybrid_evolver.chainreach(automaton,initial_set,bounding_set)

print "Exporting to postscript output...",
epsbb=Rectangle("[-8.1,8.1]x[-8.1,8.1]") # eps bounding box
eps=EpsPlot("affine_hybrid_automaton-1.eps",epsbb)
eps.set_line_style(True)
eps.set_fill_colour("cyan")
eps.write(jump_set)
eps.set_fill_colour("yellow")
eps.write(chainreach_set[mode_id])
eps.set_fill_colour("green")
eps.write(continuous_chainreach_set[mode_id])
eps.set_fill_colour("white")
eps.write(discrete_step_set[mode_id])
eps.set_fill_colour("blue")
eps.write(initial_set[mode_id])
eps.close()
print " done."