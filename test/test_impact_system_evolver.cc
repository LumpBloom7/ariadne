/***************************************************************************
 *            test_impact_system_evolver.cc
 *
 *  Copyright  2008  Pieter Collins
 *
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

#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "differentiation/sparse_differential.h"
#include "function/function_interface.h"
#include "function/affine_function.h"
#include "geometry/zonotope.h"
#include "geometry/taylor_set.h"
#include "system/impact_system.h"
#include "evaluation/impact_system_evolver.h"
#include "output/epsstream.h"


using namespace Ariadne;

int main() {
  typedef Float64 R;
  typedef ApproximateFloat64 A;

  EvolutionParameters<R> evolution_parameters;
  evolution_parameters.set_maximum_step_size(0.25);
  Evolver<ImpactSystem<R>,Zonotope<R> > evolver(evolution_parameters);
  std::cout << "evolver = " << evolver << std::endl;

  evolver.verbosity=5;

  FunctionInterface<R>* dynamic=new AffineFunction<R>(Vector<R>("[1,0]"),Matrix<R>("[0.0,0.0;-1.0,0.0]"));
  FunctionInterface<R>* reset=new AffineFunction<R>(Vector<R>("[-2,-1.0]"),Matrix<R>("[0.9,0;0,0.9]"));
  FunctionInterface<R>* guard=new AffineFunction<R>(Vector<R>("[-1]"),Matrix<R>("[0,1]"));

  ImpactSystem<R> impact_oscillator(*dynamic,*guard,*reset);

  Box<R> initial_box("[-1.03,-1.01]x[0.49,0.51]"); // initial state
  //Box<R> initial_box("[-1.11,-1.09]x[0.49,0.51]"); // initial state
  //Box<R> initial_box("[-1.21,-1.19]x[0.49,0.51]"); // initial state
  //Box<R> initial_box("[-1.51,-1.49]x[0.49,0.51]"); // initial state
  //Box<R> initial_box("[-1.81,-1.79]x[0.49,0.51]"); // initial state
  //Box<R> initial_box("[-1.83,-1.81]x[0.49,0.51]"); // initial state
  //Box<R> initial_box("[-1.86,-1.84]x[0.49,0.51]"); // initial state
  Zonotope<R> initial_set(initial_box);
  Rational time=1.5;
  Integer steps=4;
  HybridTime hybrid_time(time,steps);

  ListSet< Zonotope<R> > evolve_sets;
  ListSet< Zonotope<R> > reach_sets;
  ListSet< Zonotope<R> > intermediate_sets;

  evolver.evolution(evolve_sets,reach_sets,intermediate_sets,impact_oscillator,initial_set,hybrid_time);

  std::cout << "\nfinal:\n";
  for(uint i=0; i!=evolve_sets.size(); ++i) {
    std::cout << evolve_sets[i] << std::endl;
  }
  std::cout << "\nreach:\n";
  for(uint i=0; i!=reach_sets.size(); ++i) {
    std::cout << reach_sets[i] << std::endl;
  }
  std::cout << "\nintermediate:\n";
  for(uint i=0; i!=intermediate_sets.size(); ++i) {
    std::cout << intermediate_sets[i] << std::endl;
  }
  std::cout << "\ninitial:\n";
  std::cout << initial_set << std::endl;


  epsfstream eps;
  eps.open("test_impact_system_evolver.eps",Box<R>("[-2.625,1.625]x[-1.25,1.25]"));
  //eps << fill_colour(cyan) << Box<R>("[1,1.25]x[-2,2]");
  eps << fill_colour(cyan) << Box<R>("[-3,3]x[1,2]");
  for(uint i=0; i!=reach_sets.size(); ++i) {
    eps << fill_colour(green) << reach_sets[i];
  }
  for(uint i=0; i!=intermediate_sets.size(); ++i) {
    eps << fill_colour(yellow) << intermediate_sets[i];
  }
  for(uint i=0; i!=evolve_sets.size(); ++i) {
    eps << fill_colour(blue) << evolve_sets[i];
  }
  eps << fill_colour(blue) << initial_set;
  eps.close();

}

