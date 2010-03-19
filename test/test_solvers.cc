/***************************************************************************
 *            test_solvers.cc
 *
 *  Copyright  2008-10  Pieter Collins
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

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "solver.h"
#include "function.h"
#include "taylor_function.h"
#include "vector.h"

#include "test.h"

using namespace Ariadne;
using namespace std;

typedef Vector<Interval> IntervalVector;

class TestSolver
{
  private:
    scoped_ptr<SolverInterface> solver;
  public:
    TestSolver(const SolverInterface& s)
        : solver(s.clone()) { }

    int test() {
        ARIADNE_TEST_CALL(test_solve());
        ARIADNE_TEST_CALL(test_implicit());
        return 0;
    }

    void test_solve() {
        ScalarFunction x=ScalarFunction::coordinate(1,0);
        IntervalVector d(1, Interval(0.0,1.0));
        VectorFunction f(1u, (x*x+1)*x-1);
        IntervalVector p=solver->solve(f,d);
        ARIADNE_TEST_BINARY_PREDICATE(subset,p[0],Interval(0.6823,0.6824));
    }

    void test_implicit() {
        TaylorModel::set_default_sweep_threshold(1e-12);

        ScalarFunction aa=ScalarFunction::coordinate(1,0);
        ScalarFunction a=ScalarFunction::coordinate(2,0);
        ScalarFunction x=ScalarFunction::coordinate(2,1);
        ScalarFunction bb;
        IntervalVector p,r;
        VectorFunction f;
        VectorTaylorFunction h,e;

        // Test solution of x-a=0. This should be very easy to solve.
        p=IntervalVector(1, Interval(-0.25,0.25));
        r=IntervalVector(1, Interval(-2.0,2.0));
        f=VectorFunction(1u,x-a);
        h=solver->implicit(f,p,r);
        e=VectorTaylorFunction(p,VectorFunction(1u,aa));
        ARIADNE_TEST_COMPARE(norm((h-e).range()),<,1e-8);

        // Test solution of 4x^2+x-4-a=0 on [0.875,1.125]. There is a unique solution with positive derivative.
        p=IntervalVector(1, Interval(0.875,1.125));
        r=IntervalVector(1, Interval(0.25,1.25));
        f=VectorFunction(1u,(x*x+1)*x-a);
        h=solver->implicit(f,p,r);
        bb=ScalarFunction(aa-p[0].midpoint())/p[0].radius();
        e=VectorTaylorFunction( p, VectorFunction( 1u, 0.682328+bb*(0.0521547+bb*(-0.0023232+bb*0.000147778)) ) );
        ARIADNE_TEST_COMPARE(norm((h-e).range()),<,1e-4);

        // Test solution of 4x^2+x-4-a=0 on [-0.25,0.25]. There is a unique solution with positive derivative.
        p=IntervalVector(1, Interval(-0.25,0.25));
        r=IntervalVector(1, Interval(0.25,2.0));
        f=VectorFunction(1u,4*x+x*x-a-4);
        h=solver->implicit(f,p,r);
        bb=ScalarFunction(aa-p[0].midpoint())/p[0].radius();
        e=VectorTaylorFunction( p, VectorFunction( 1u, 0.828427+bb*(0.0441942+bb*(-0.000345267+bb*0.00000539468)) ) );
        ARIADNE_TEST_COMPARE(norm((h-e).range()),<,1e-4);
    }

};


int main() {
    KrawczykSolver krawczyk_solver(1e-5,12);
    krawczyk_solver.verbosity=0;
    TestSolver(krawczyk_solver).test();

    FactoredKrawczykSolver factored_krawczyk_solver(1e-5,12);
    factored_krawczyk_solver.verbosity=0;
    TestSolver(factored_krawczyk_solver).test();
    std::cerr<<"INCOMPLETE "<<std::flush;

}
