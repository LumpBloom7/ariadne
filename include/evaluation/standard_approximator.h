/***************************************************************************
 *            standard_approximator.h
 *
 *  Copyright  2007  Pieter Collins
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
 
/*! \file standard_approximator.h
 *  \brief Methods for approximating basic sets
 */

#ifndef ARIADNE_STANDARD_APPROXIMATOR_H
#define ARIADNE_STANDARD_APPROXIMATOR_H

#include "base/types.h"
#include "base/declarations.h"
#include "numeric/declarations.h"
#include "linear_algebra/declarations.h"
#include "geometry/declarations.h"

#include "approximator_interface.h"

namespace Ariadne {
  namespace Evaluation {

    /*! \brief Geomerical approximation schemes.
     *  \ingroup Approximation
     */
    template<class BS>
    class StandardApproximator
      : public ApproximatorInterface<BS>
    {
      typedef typename BS::real_type R;
      typedef Numeric::Interval<R> I;
     public:
      //virtual ~StandardApproximator();
      StandardApproximator();
      StandardApproximator(const StandardApproximator<BS>& approx);
      virtual StandardApproximator<BS>* clone() const;
      virtual BS over_approximation(const Geometry::Box<R>& r) const;
      virtual Geometry::Box<R> bounding_box(const BS& bs) const;
      virtual Geometry::GridCellListSet<R> outer_approximation(const BS& bs, const Geometry::Grid<R>& g) const;
      virtual std::pair<BS,BS> subdivide(const BS&) const;
    };



  }
}

#endif /* ARIADNE_STANDARD_APPROXIMATOR_H */
