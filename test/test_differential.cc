/***************************************************************************
 *            test_differential.cc
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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

#include <cassert>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>

#include "test_float.h"
#include "numeric/rational.h"
#include "numeric/differential.h"
#include "linear_algebra/vector.h"

#include "test.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace std;

template class Differential< Float, Vector<Float> >;
template class Differential< Rational, Vector<Rational> >;

template<class X, class V>
class TestDifferential {
 private:
  Differential<X,V> x1,x2,x3,x4,x5;
 public:
  TestDifferential() {
    x1=Differential<X,V>(2.0,2,0);
    x2=Differential<X,V>(2.0,2,1);
    x3=Differential<X,V>(2.0,"[3,2]");
    x4=Differential<X,V>(-1.0,"[1,4]");
    x5=Differential<X,V>(2,"[1,-1]");

    test_add();
    test_sub();
    test_mul();
    test_div();
    test_pow();
  }

  void test_add() {
    cout << x1 << "+" << x2 << " = " << x1+x2 << std::endl;
    assert(((x1+x2)==Differential<X,V>("4","[1,1]")));
  }

  void test_sub() {
    cout << x1 << "-" << x2 << " = " << x1-x2 << std::endl;
    assert(((x1-x2)==Differential<X,V>("0","[1,-1]")));
  }

  void test_mul() {
    cout << x3 << "*" << x4 << " = " << x3*x4 << std::endl;
    assert(((x3*x4)==Differential<X,V>("-2","[-1,6]")));
  }

  void test_div() {
    cout << x3 << "/" << x4 << " = " << x3/x4 << std::endl;
    assert(((x3/x4)==Differential<X,V>("-2","[-5,-10]")));
  }

  void test_pow() {
    cout << x5 << "^5 = " << pow(x5,5) << std::endl;
    assert((pow(x5,5)==Differential<X,V>("32","[80,-80]")));
  }

};

template<class X>
class TestScalarDerivative {
 private:
  ScalarDerivative<X> x1,x2,x3;
 public:
  TestScalarDerivative() {
    x1=ScalarDerivative<X>(3,1.0,1);
    x2=ScalarDerivative<X>(3,2.0,1);
    x3=ScalarDerivative<X>("[2.0,1.0,0.1,0.0]");
    
    test_degree();
    test_add();
    test_sub();
    test_mul();
    test_div();
    test_pow();
  }

  void test_degree() {
    assert(x1.degree()==3);
  }

  void test_add() {
    cout << x1 << "+" << x2 << " = " << x1+x2 << std::endl;
    assert((x1+x2)==ScalarDerivative<X>("[3,2,0,0]"));
  }

  void test_sub() {
    cout << x1 << "-" << x2 << " = " << x1-x2 << std::endl;
    assert((x1-x2)==ScalarDerivative<X>("[-1,0,0,0]"));
  }

  void test_mul() {
    cout << x1 << "*" << x2 << " = " << x1*x2 << std::endl;
    assert((x1*x2)==ScalarDerivative<X>("[2,3,2,0]"));
  }

  void test_div() {
    ScalarDerivative<X> x3("[2,3,4]");
    ScalarDerivative<X> x4("[1,0,0]");
    cout << x3 << "/" << x4 << " = " << x3/x4 << std::endl;
    assert((x3/x4)==x3);
    cout << x4 << "/" << x3 << " = " << x4/x3 << std::endl;
    assert((x4/x3)==ScalarDerivative<X>("[0.5,-0.75,1.25]"));
    cout << 1 << "/" << x2 << " = " << 1/x2 << std::endl;
    assert((1/x2)==ScalarDerivative<X>("[0.5,-0.25,0.25,-0.375]"));
    cout << x1 << "/" << x2 << " = " << x1/x2 << std::endl;
    assert((x1/x2)==ScalarDerivative<X>("[0.5,0.25,-0.25,0.375]"));
  }

  void test_pow() {
    cout << x2 << "^5 = " << pow(x2,5) << std::endl;
    assert(pow(x2,5)==ScalarDerivative<X>("[32,80,160,240]"));
  }

};


int main() {
  TestDifferential<Rational, Vector<Rational> > t1;
  cout << endl;
  TestScalarDerivative<Rational> t2;
  
  return 0;
}