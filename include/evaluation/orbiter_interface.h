/***************************************************************************
 *            orbiter_interface.h
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
 
/*! \file orbiter_interface.h
 *  \brief Methods for computing orbits of boxes.
 */

#ifndef ARIADNE_ORBITER_INTERFACE_H
#define ARIADNE_ORBITER_INTERFACE_H

#include <boost/smart_ptr.hpp>

#include "base/types.h"
#include "base/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"


namespace Ariadne {
  namespace Evaluation {


    /*! \brief A class for computing the evolution of a continuous-time autonomous system.
     *  \ingroup Integrators
     */
    template<class R>
    class VectorFieldOrbiterInterface 
    {
     public:
      /*! \brief Destructor. */
     virtual ~VectorFieldOrbiterInterface();

      /*! \brief Make a dynamically-allocated copy. */
      virtual VectorFieldOrbiterInterface<R>* clone() const = 0;

      /*! \brief Compute the orbit of a box under time \a t of a vector field. */
      virtual 
      Geometry::OrbitInterface<Numeric::Rational>*
      orbit(const System::VectorFieldInterface<R>& vf, const Geometry::Box<R>& bx, const Numeric::Rational& t) const = 0;

       
    };

    template<class R> VectorFieldOrbiterInterface<R>::~VectorFieldOrbiterInterface() { }


  }

}

#endif /* ARIADNE_ORBITER_INTERFACE_H */
