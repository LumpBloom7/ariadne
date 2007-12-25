/***************************************************************************
 *            test_constraint_based_hybrid_evolver.cc
 *
 *  Copyright  2007  Pieter Collins
 *  Email  Pieter.Collins@cwi.nl
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
#include "test_float.h"
#include "geometry/set_interface.h"
#include "geometry/set_reference.h"
#include "geometry/hybrid_set.h"
#include "geometry/rectangle.h"
#include "geometry/empty_set.h"
#include "geometry/polyhedral_set.h"
#include "geometry/linear_constraint.h"
#include "system/affine_map.h"
#include "system/affine_vector_field.h"
#include "system/constraint_based_hybrid_automaton.h"
#include "evaluation/evolution_parameters.h"
#include "evaluation/applicator.h"
#include "evaluation/lohner_integrator.h"
#include "evaluation/affine_integrator.h"
#include "evaluation/constraint_based_hybrid_evolver.h"
#include "evaluation/constraint_based_hybrid_evolver.h"
#include "output/epsstream.h"
#include "output/logging.h"

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

static const DiscreteState dstate1(1);
static const DiscreteState dstate2(2);
static const DiscreteState dstate3(3);
static const DiscreteEvent event20_id(6);
static const DiscreteEvent event12_id(4);
static const DiscreteEvent event23_id(5);
  

template<class R>
ConstraintBasedHybridEvolver<R> 
construct_evolver() 
{
  typedef Interval<R> I;
  typedef Zonotope<I,I> BS;

  EvolutionParameters<R> parameters;
  parameters.set_maximum_step_size(0.125);
  parameters.set_lock_to_grid_time(0.5);
  parameters.set_maximum_basic_set_radius(0.25);
  parameters.set_grid_length(0.125);
  
  Applicator<R> applicator;
  C1LohnerIntegrator<R> lohner_integrator; 
  const ApplicatorInterface<BS>& applicator_interface=applicator;
  const IntegratorInterface<BS>& integrator_interface=lohner_integrator;
  return ConstraintBasedHybridEvolver<R>(parameters,applicator_interface,integrator_interface);
}



template<class R>
ConstraintBasedHybridAutomaton<R> 
construct_automaton() 
{
    AffineVectorField<R> dynamic1(Matrix<R>("[-1,-0.25;0.25,-1]"),Vector<R>("[5,0]"));
    AffineVectorField<R> dynamic2(Matrix<R>("[0,0; 0,0]"),Vector<R>("[0.75,1]"));
    AffineVectorField<R> dynamic3(Matrix<R>("[0,0; 0,0]"),Vector<R>("[1,-0.75]"));

    AffineMap<R> reset12(Matrix<R>("[0,1;-1,0]"),Vector<R>("[-2,3]"));
    LinearConstraint<R> guard12(Vector<R>("[1,0]"),Geometry::less,R(2));
    AffineMap<R> reset23(Matrix<R>("[1,0;0,1]"),Vector<R>("[0,-3]"));
    LinearConstraint<R> activation23(Vector<R>("[0,1]"),Geometry::less,R(2));
    LinearConstraint<R> invariant2(Vector<R>("[0,1]"),Geometry::less,R(2.5));
    
    ConstraintBasedHybridAutomaton<R> automaton("");

    automaton.new_mode(dstate1,dynamic1);
    automaton.new_mode(dstate2,dynamic2);
    automaton.new_mode(dstate3,dynamic3);
    automaton.new_invariant(event20_id,dstate2,invariant2);
    automaton.new_forced_transition(event12_id,dstate1,dstate2,reset12,guard12);
    automaton.new_unforced_transition(event23_id,dstate2,dstate3,reset23,activation23);

    cout << "automaton = " << flush;
    cout << automaton << endl << endl;

    return automaton;
}





template<class R>
class TestConstraintBasedHybridEvolver
{
 public:
  typedef Interval<R> I;
  typedef Zonotope<I> BS;

  ConstraintBasedHybridAutomaton<R> automaton;
  ConstraintBasedHybridEvolver<R> evolver;
  Polyhedron<R> guard;
  Polyhedron<R> activation;
  Polyhedron<R> invariant;
    
  TestConstraintBasedHybridEvolver()
    : automaton(construct_automaton<R>()),
      evolver(construct_evolver<R>())
  {
  }


  int test_evolution() {
    time_type tr=3.0;
    time_type t=1.625;
    //time_type t=1.75;
    t=1.90;
    size_type n=2;

    DiscreteState initial_discrete_mode = dstate1;
    Zonotope<I,I> initial_basic_set(Point<R>("(0.5,0)"),Matrix<R>("[0.03125,0.0; 0.0,0.03125]"));
    HybridListSet< Zonotope<I,I> > initial_set(automaton.locations());
    initial_set.adjoin(initial_discrete_mode,initial_basic_set);
  

    
    cout << endl;

    Box<R> bounding_box("[-3,3]x[-3,3]");
    Polyhedron<R> bounding_polyhedron(bounding_box);
    Polyhedron<R> guard_polyhedron("[[-1,0;-2]]");
    Polyhedron<R> activation_polyhedron("[[0,-1;-2]]");
    
    HybridListSet< Zonotope<I,I> > evolved_set=initial_set;
    evolved_set.clear();
    HybridListSet< Zonotope<I,I> > reached_set=initial_set;
    reached_set.clear();
    try {
      evolved_set=evolver.upper_evolve(automaton,initial_set,t,n);
      reached_set=evolver.upper_reach(automaton,initial_set,tr,n);
    } 
    catch(const std::exception& e) {
      cout << "\nCaught: " << e.what() << endl;
    }
    cout << "\n" << endl;
    cout << "trace=" << evolver.trace() << endl << endl;

    cout << "initial_set = " << initial_set << endl;
    cout << "evolved_set = " << evolved_set << endl;
    cout << "reached_set = " << reached_set << endl;
    
    epsfstream eps;
    eps.open("test_hybrid_evolution-trace.eps",bounding_box);
    eps << fill_colour(magenta) << closed_intersection(bounding_polyhedron,guard_polyhedron);
    eps << fill_colour(cyan) << closed_intersection(bounding_polyhedron,activation_polyhedron);

    for(uint i=0; i!=evolver.trace().size(); ++i) {
      switch(evolver.trace()[i].discrete_state().id()) {
        case 1: eps << fill_colour(Output::green); break;
        case 2: eps << fill_colour(Output::red); break;
        case 3: eps << fill_colour(Output::cyan); break;
      }
      eps << evolver.trace()[i].continuous_state_set();
    }      
    eps.close();

    eps.open("test_hybrid_evolution-evolve.eps",bounding_box);
    eps << fill_colour(magenta) << closed_intersection(bounding_polyhedron,guard_polyhedron);
    eps << fill_colour(cyan) << closed_intersection(bounding_polyhedron,activation_polyhedron);
    eps << fill_colour(Output::yellow) << initial_set[dstate1];
    eps << fill_colour(Output::white) << evolved_set[dstate1];
    eps << fill_colour(Output::white) << evolved_set[dstate2];
    eps << fill_colour(Output::white) << evolved_set[dstate3];
    eps.close();
    
    eps.open("test_hybrid_evolution-reach.eps",bounding_box);
    eps << fill_colour(magenta) << closed_intersection(bounding_polyhedron,guard_polyhedron);
    eps << fill_colour(cyan) << closed_intersection(bounding_polyhedron,activation_polyhedron);
    eps << fill_colour(Output::green) << reached_set[dstate1];
    eps << fill_colour(Output::red) << reached_set[dstate2];
    eps << fill_colour(Output::blue) << reached_set[dstate3];
    eps << fill_colour(Output::yellow) << initial_set[dstate1];
    eps.close();

    return 0;
  }
  

  int test() {
    test_evolution();
    return 0;
  }

};



int main(int nargs, const char* args[]) {
  int verbosity = 0;
  if(nargs>1) { verbosity=std::atoi(args[1]); }
  set_hybrid_evolver_verbosity(verbosity);

  TestConstraintBasedHybridEvolver<Flt>().test();
  cerr << "INCOMPLETE ";
  return 0;
}
