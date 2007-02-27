/***************************************************************************
 *            python/export_models.cc
 *
 *  Copyright  2005-6  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is diself_ns::stributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "real_typedef.h"

#include "models/henon.h"
#include "models/duffing.h"
#include "models/vanderpol.h"

#include "models/lorenz.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;

#include <boost/python.hpp>
using namespace boost::python;
  

template<class R>
void export_henon_map() {
  typedef typename Numeric::traits<R>::arithmetic_type F;
  
  class_< HenonMap<R>, bases< Map<R> > >("HenonMap",init<R,R>())
    .def("argument_dimension", &HenonMap<R>::argument_dimension)
    .def("result_dimension", &HenonMap<R>::result_dimension)
    .def("__call__", (Point<F>(HenonMap<R>::*)(const Point<F>&)const)(&HenonMap<R>::image))
    .def("jacobian", (Matrix<F>(HenonMap<R>::*)(const Point<F>&)const)(&HenonMap<R>::jacobian))
    .def(self_ns::str(self))
  ;
}


template<class R>
void export_duffing_equation() {
  typedef typename Numeric::traits<R>::arithmetic_type F;
  
  class_< DuffingEquation<R>, bases< VectorField<R> > >("DuffingEquation",init<R,R>())
    .def("dimension", &DuffingEquation<R>::dimension)
    .def("__call__", (Vector<F>(DuffingEquation<R>::*)(const Point<F>&)const)(&DuffingEquation<R>::image))
    .def("jacobian", (Matrix<F>(DuffingEquation<R>::*)(const Point<F>&)const)(&DuffingEquation<R>::jacobian))
    .def(self_ns::str(self))
  ;
}

template<class R>
void export_van_der_pol_equation() {
  typedef typename Numeric::traits<R>::arithmetic_type F;
  
  class_< VanDerPolEquation<R>, bases< VectorField<R> > >("VanDerPolEquation",init<R>())
    .def("dimension", &VanDerPolEquation<R>::dimension)
    .def("__call__", (Vector<F>(VanDerPolEquation<R>::*)(const Point<F>&)const)(&VanDerPolEquation<R>::image))
    .def("jacobian", (Matrix<F>(VanDerPolEquation<R>::*)(const Point<F>&)const)(&VanDerPolEquation<R>::jacobian))
    .def(self_ns::str(self))
  ;
}


template<class R>
void export_lorenz_system() 
{
  typedef typename Numeric::traits<R>::arithmetic_type F;
  
  class_< LorenzSystem<R>, bases< VectorField<R> > >("LorenzSystem",init<R,R,R>())
    .def("dimension", &LorenzSystem<R>::dimension)
    .def("smoothness", &LorenzSystem<R>::smoothness)
    .def("__call__", (Vector<F>(LorenzSystem<R>::*)(const Point<F>&)const)(&LorenzSystem<R>::image))
    .def("jacobian", (Matrix<F>(LorenzSystem<R>::*)(const Point<F>&)const)(&LorenzSystem<R>::jacobian))
    .def(self_ns::str(self))
  ;
}

template void export_henon_map<Real>();
template void export_duffing_equation<Real>();
template void export_van_der_pol_equation<Real>();
template void export_lorenz_system<Real>();