/***************************************************************************
 *            python/export_Matrix.cc
 *
 *  17 November 2005
 *  Copyright  2005  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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

#include <utility>  //for std::pair


#include "numeric/interval.h"
#include "numeric/integer.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "linear_algebra/interval_matrix.h"

#include "python/typedefs.h"
#include "python/python_utilities.h"
using namespace Ariadne;

#include <boost/python.hpp>
using namespace boost::python;

inline Field fmatrix_getitem(const FMatrix& A, tuple index) {
  uint i=extract<uint>(index[0]);
  uint j=extract<uint>(index[1]);
  return A(i,j);
}

inline void fmatrix_setitem(FMatrix& A, tuple index, Field x) {
  uint i=extract<uint>(index[0]);
  uint j=extract<uint>(index[1]);
  A(i,j)=x;
}

inline void fmatrix_setitem_from_real(FMatrix& A, tuple index, Real x) {
  uint i=extract<uint>(index[0]);
  uint j=extract<uint>(index[1]);
  A(i,j)=Ariadne::convert_to<Field>(x);
}

inline void fmatrix_setitem_from_double(FMatrix& A, tuple index, double x) {
  uint i=extract<uint>(index[0]);
  uint j=extract<uint>(index[1]);
  A(i,j)=Ariadne::convert_to<Field>(x);
}

inline Real rmatrix_getitem(const RMatrix& A, tuple index) {
  uint i=extract<uint>(index[0]);
  uint j=extract<uint>(index[1]);
  return A(i,j);
}

inline void rmatrix_setitem(RMatrix& A, tuple index, Real x) {
  uint i=extract<uint>(index[0]);
  uint j=extract<uint>(index[1]);
  A(i,j)=x;
}

inline void rmatrix_setitem_from_double(RMatrix& A, tuple index, double x) {
  uint i=extract<uint>(index[0]);
  uint j=extract<uint>(index[1]);
  A(i,j)=Ariadne::convert_to<Real>(x);
}

inline RVector rmatrix_vprod(const RMatrix& A, const RVector v) {
  return prod(A,v);
}

inline RInterval iMatrix_getitem(const RIntervalMatrix& A, tuple index) {
  uint i=extract<uint>(index[0]);
  uint j=extract<uint>(index[1]);
  return A(i,j);
}

inline void iMatrix_setitem(RIntervalMatrix& A, tuple index, RInterval x) {
  uint i=extract<uint>(index[0]);
  uint j=extract<uint>(index[1]);
  A(i,j)=x;
}

inline FMatrix rmatrix_inverse(const RMatrix& A) {
  return LinearAlgebra::inverse(A);
}

void export_Matrix() {
  class_<RMatrix>("Matrix",init<int,int>())
    .def(init<std::string>())
    .def(init<RMatrix>())
    .def("__mul__",&rmatrix_vprod)
    .def("__getitem__",&rmatrix_getitem)
    .def("__setitem__",&rmatrix_setitem)
    .def("__setitem__",&rmatrix_setitem_from_double)
    .def(self_ns::str(self))    // __self_ns::str__
  ;
    
  def("inverse",&rmatrix_inverse);
  def("exp_approx",&LinearAlgebra::exp_approx<Real>);

#ifndef REAL_IS_A_FIELD 
  class_<FMatrix>("RationalMatrix",init<int,int>())
    .def(init<FMatrix>())
    .def("__getitem__",&fmatrix_getitem)
    .def("__setitem__",&fmatrix_setitem)
    .def("__setitem__",&fmatrix_setitem_from_real)
    .def("__setitem__",&fmatrix_setitem_from_double)
    .def(self_ns::str(self))    // __self_ns::str__
  ;
#endif
}

void export_IntervalMatrix() {
  class_<RIntervalMatrix>("IntervalMatrix",init<int,int>())
    .def(init<RIntervalMatrix>())
    .def("__getitem__",&iMatrix_getitem)
    .def("__setitem__",&iMatrix_setitem)
    .def(self_ns::str(self))    // __self_ns::str__
  ;
}
