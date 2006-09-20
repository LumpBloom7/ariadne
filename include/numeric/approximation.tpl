/***************************************************************************
 *            approximation.tpl
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
 
/*! \file approximation.h
 *  \brief Approximation of numerical types.
 */

#ifndef _ARIADNE_APPROXIMATION_H
#define _ARIADNE_APPROXIMATION_H

#include "approximation.h"

#include "float64.h"
#include "mpfloat.h"
#include "dyadic.h"
#include "rational.h"

namespace Ariadne {
  namespace Numeric {
    
    /*! \brief Approximate \a x by an element of \p Res. */
    template<typename Res, typename Arg> Res approximate_by(const Arg& x);
        
    template<> double approximate_by(const double& x) {
      return x;
    }
    template<> double approximate_by(const MPFloat& x) { 
      return x.get_d(); 
    }
    
    template<> double approximate_by(const Dyadic& x) {
      return x.get_d();
    }
    
    template<> double approximate_by(const Rational& q) {
      return q.get_d();
    }
    
    template<> MPFloat approximate_by(const Rational& q) {
      return MPFloat(q.get_d());
    }
        
    template<> double approximate(const double& x, const double& e) {
      return x;
    }
    
    template<> Rational approximate(const Rational& q, const Rational& e) {
      return q;
    }
    
    template<> Dyadic approximate(const Rational& q, const Rational& e) {
      Dyadic x=Dyadic(q);
      assert(abs(Rational(x)-q)<e);
      return x;
    }
      
    template<> Dyadic approximate(const Rational& q, const Dyadic& e) {
      return approximate<Dyadic>(q,Rational(e));
    }
      
    template<> MPFloat approximate(const Rational& q, const Rational& e) {
      MPFloat x=MPFloat(mpf_class(q));
      assert(abs(Rational(mpf_class(x))-q)<e);
      return x;
    }
      
    template<> MPFloat approximate(const Rational& q, const MPFloat& e) {
      return approximate<MPFloat>(q,Rational(mpf_class(e)));
    }
  }    
}

#endif /* _ARIADNE_APPROXIMATION_H */