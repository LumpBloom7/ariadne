/***************************************************************************
 *            test_hybrid_evolution.cc
 *
 *  Copyright  2006  Alberto Casagrande,  Pieter Collins
 *  Email  casagrande@dimi.uniud.it  Pieter.Collins@cwi.nl
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include <iostream>
#include <fstream>
#include <string>

#include "ariadne.h"
#include "debug.h"
#include "real_typedef.h"
#include "geometry/set.h"
#include "geometry/hybrid_set.h"
#include "geometry/rectangle.h"
#include "geometry/polyhedral_set.h"
#include "system/affine_map.h"
#include "system/affine_vector_field.h"
#include "system/hybrid_automaton.h"
#include "evaluation/applicator.h"
#include "evaluation/lohner_integrator.h"
#include "evaluation/affine_integrator.h"
#include "evaluation/hybrid_evolver.h"
#include "output/epsfstream.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Combinatoric;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace Ariadne::Evaluation;
using namespace Ariadne::Output;
using namespace std;

template<class R> int test_hybrid_evolution();
  
namespace Ariadne { namespace Evaluation { extern int verbosity; } }

int main() {
  test_hybrid_evolution<Real>();
  cerr << "INCOMPLETE ";
  return 0;
}

template<class R>
int test_hybrid_evolution() 
{  
  Evaluation::verbosity=7; 
  Rectangle<R> r("[-7.5,7.5]x[-7.5,7.5]");
  Polyhedron<R> p(r);

  AffineVectorField<R> dynamic(Matrix<R>("[-2,-1;1,-2]"),Vector<R>("[-1,0]"));
  AffineMap<R> reset(Matrix<R>("[1,0;0,-1]"),Vector<R>("[0,0]"));
  
  PolyhedralSet<R> invariant(r);
  PolyhedralSet<R> activation(Rectangle<R>("[-7.5,7.5]x[-3,-2]"));
  
  HybridAutomaton<R> automaton("Affine automaton");
  id_type mode_id=2;
  const DiscreteMode<R>& mode=automaton.new_mode(mode_id,dynamic,invariant);
  id_type event_id=3;
  const DiscreteTransition<R>& transition=automaton.new_transition(event_id,mode_id,mode_id,reset,activation);
  
  cout << mode << endl << transition << endl;
  
  Applicator<R> apply;
  //LohnerIntegrator<R> integrator(0.1,0.1,0.1);
  AffineIntegrator<R> integrator(0.125,0.5,0.25); 
  HybridEvolver<R> hybrid_evolver(apply,integrator);
  
  RegularGrid<R> grid(Vector<R>("[0.25,0.25]"));
  FiniteGrid<R> finite_grid(grid,LatticeBlock("[-32,32]x[-32,32]"));
  Rectangle<R> bounding_box=finite_grid.extent();
  
  Rectangle<R> initial_rectangle("[-6.96875,-6.9375]x[-6.96875,-6.9375]");
  HybridGridMaskSet<R> initial_set;
  initial_set.new_location(mode_id,finite_grid);
  initial_set[mode_id].adjoin_over_approximation(initial_rectangle);
  cout << "initial_set.discrete_locations()=" << initial_set.locations() << endl;

  HybridGridMaskSet<R> bounding_set;
  bounding_set.new_location(mode_id,finite_grid);
  bounding_set[mode_id].adjoin_over_approximation(bounding_box);
  cout << "bounding_set.discrete_locations()=" << bounding_set.locations() << endl;

  HybridGridMaskSet<R> continuous_chainreach=hybrid_evolver.continuous_chainreach(automaton,initial_set,bounding_set);
  cout << "Reached " << continuous_chainreach[mode_id].size() << " cells out of " << continuous_chainreach[mode_id].capacity() 
       << " by continuous evolution" << endl << endl;
  
  epsfstream eps;
  eps.open("test_hybrid_evolution-1.eps",bounding_box.expand(0.5));
  eps.set_fill_colour("white");
  eps << bounding_box;
  eps.set_fill_colour("green");
  eps << continuous_chainreach[mode_id];
  eps.set_fill_colour("blue");
  eps << initial_set[mode_id];
  eps.close();


  HybridGridCellListSet<R> initial_activated_set;
  initial_activated_set.new_location(mode_id,grid);
  initial_activated_set[mode_id].adjoin(over_approximation(Rectangle<R>("[-6,-5]x[-3.5,-1.0]"),grid));
  cout << "initial_activated_set[0].size()=" << initial_activated_set[mode_id].size() << " cells" << endl;
  HybridGridCellListSet<R> discrete_reach=hybrid_evolver.discrete_step(automaton,initial_activated_set);
  cout << "Reached " << discrete_reach[mode_id].size() << " cells by discrete step" << endl << endl;
  
  eps.open("test_hybrid_evolution-2.eps",bounding_box.expand(0.5));
  eps.set_fill_colour("white");
  eps << bounding_box;
  eps.set_fill_colour("cyan");
  eps << static_cast<Polyhedron<R>&>(activation);
  eps.set_fill_colour("green");
  eps << discrete_reach[mode_id];
  eps.set_fill_colour("blue");
  eps << initial_activated_set[mode_id];
  eps.close();
  
  
  cout << "bounding_set[0].size()=" << bounding_set[mode_id].size() << " cells out of " << bounding_set[mode_id].capacity() << endl;
  assert(subset(initial_set[mode_id],bounding_set[mode_id]));

  cout << "Here" << endl;
  HybridGridMaskSet<R> chainreach=hybrid_evolver.chainreach(automaton,initial_set,bounding_set);
  cout << "Here" << endl;
  cout << "Reached " << chainreach[mode_id].size() << " cells out of " << chainreach[mode_id].capacity() << endl;
  
  eps.open("test_hybrid_evolution-3.eps",bounding_box.expand(0.5));
  eps.set_fill_colour("white");
  eps << bounding_box;
  eps.set_fill_colour("cyan");
  eps << static_cast<Polyhedron<R>&>(activation);
  eps.set_line_style(true);
  eps.set_fill_colour("yellow");
  eps << chainreach[mode_id];
  eps.set_fill_colour("blue");
  eps << initial_set[mode_id];
  eps.close();
  
  
  return 0;
}