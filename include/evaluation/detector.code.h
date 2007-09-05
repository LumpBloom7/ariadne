/***************************************************************************
 *            detector.code.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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
 
#include "system/vector_field_interface.h"
#include "geometry/constraint_interface.h"
#include "geometry/basic_set_interface.h"

namespace Ariadne {

template<class R>
Evaluation::Detector<R>::~Detector()
{
}

template<class R>
Evaluation::Detector<R>::Detector(const Detector<R>& det) 
{
}

template<class R>
Evaluation::Detector<R>*
Evaluation::Detector<R>::clone() const
{
  return new Detector<R>(*this);
}

template<class R>
Numeric::Interval<R> 
Evaluation::Detector<R>::value(const Geometry::ConstraintInterface<R>& c, 
                               const Geometry::BasicSetInterface<R>& bs)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}


template<class R>
Evaluation::TimeModel<R> 
Evaluation::Detector<R>::crossing_time(const System::VectorFieldInterface<R> vf, 
                                       const Geometry::ConstraintInterface<R>& c, 
                                       const Geometry::BasicSetInterface<R>& bs)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}



}