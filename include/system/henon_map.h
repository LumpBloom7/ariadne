/***************************************************************************
 *            henon_map.h
 *
 *  Wed Feb  2 18:52:36 2005
 *  Copyright  2005  Alberto Casagrande
 *  casagrande@dimi.uniud.it
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
 
/*! \file henon_map.h
 *  \brief The Henon map \f$(x,y) \rightarrow (a-x^2-by,x)\f$.
 */

#ifndef _ARIADNE_HENON_MAP_H
#define _ARIADNE_HENON_MAP_H

#include "../linear_algebra/matrix.h"
#include "../linear_algebra/interval_matrix.h"

#include "../geometry/point.h"
#include "../geometry/rectangle.h"
#include "../geometry/parallelotope.h"

#include "../system/map.h"

#include "../evaluation/apply.h"

namespace Ariadne {
  namespace System {

    /*! \brief The Henon map. */
    template <typename R>
    class HenonMap : public Map<R> 
    {
     public:
      typedef R real_type;
      typedef Geometry::Point<real_type> state_type;
      
      inline explicit HenonMap(real_type a=real_type(1.5),real_type b=real_type(0.3)) : _a(a), _b(b) { }
      
      /*! \brief  The map applied to a state. */
      virtual state_type operator() (const state_type& x) const;
      /*! \brief  The map applied to a rectangle basic set. */
      virtual Geometry::Rectangle<R> operator() (const Geometry::Rectangle<R>& r) const;
      /*! \brief  The map applied to a parallelotope basic set. */
      virtual Geometry::Parallelotope<R> operator() (const Geometry::Parallelotope<R>& r) const;
      
      /*! \brief  The derivative of the map at a point. */
      virtual LinearAlgebra::Matrix<R> derivative(const state_type& x) const;
      /*! \brief  The derivative of the map over a rectangular basic set. */
      virtual LinearAlgebra::IntervalMatrix<R> derivative(const Geometry::Rectangle<R>& r) const;
            
      /*! \brief  The parameter a. */
      const real_type& a() const { return _a; }
      /*! \brief  The parameter b. */
      const real_type& b() const { return _b; }
      
      /*! \brief  The dimension of the argument. */
      dimension_type argument_dimension() const { return 2; }
      /*! \brief The dimension of the result. */
      dimension_type result_dimension() const { return 2; }
      
      bool invertible() const {
        return _b!=0; }
        
      std::string name() const { return "HenonMap"; }
     private:
      real_type _a;
      real_type _b;
    };
      
    template <typename R>
    Geometry::Point<R>
    HenonMap<R>::operator() (const Geometry::Point<R>& p) const
    {
      Geometry::Point<R> result(2); 
      const R& x=p[0];
      const R& y=p[0];
      result[0]=_a-x*x-_b*y; 
      result[1]=x; 
      return result;
    }
     
    template <typename R>
    Geometry::Rectangle<R>
    HenonMap<R>::operator() (const Geometry::Rectangle<R>& A) const
    {
      Geometry::Rectangle<R> result(2); 
      result[0]=_a-A[0]*A[0]-_b*A[1]; 
      result[1]=A[0]; 
      return result;
    }
     
    template <typename R>
    Geometry::Parallelotope<R>
    HenonMap<R>::operator() (const Geometry::Parallelotope<R>& A) const
    {
      return Evaluation::C1Applicator<R>().apply(*this,A);
    }
     
    template <typename R>
    LinearAlgebra::Matrix<R>
    HenonMap<R>::derivative(const Geometry::Point<R>& x) const
    {
      LinearAlgebra::Matrix<R> result(2,2); 
      result(0,0) = -2*x[0];
      result(0,1) = -_b;
      result(1,0) = 1;
      result(1,1) = 0;
      return result;
    }
     
    template <typename R>
    LinearAlgebra::IntervalMatrix<R>
    HenonMap<R>::derivative(const Geometry::Rectangle<R>& r) const
    {
      LinearAlgebra::IntervalMatrix<R> result(2,2); 
      result(0,0) = R(-2)*r[0];
      result(0,1) = Interval<R>(R(-_b));
      result(1,0) = Interval<R>(R(1));
      result(1,1) = Interval<R>(R(0));
      return result;
    }
     
     
    template <typename R>
    std::ostream& operator<<(std::ostream& os, const HenonMap<R>& hm) {
      os << "HenonMap( a=" << hm.a() << ", b=" << hm.b() << " )";
      return os;
    }
    
    
    
  }
}


#endif /* _ARIADNE_HENON_MAP_H */