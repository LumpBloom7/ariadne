/***************************************************************************
 *            rational.h
 *
 *  Copyright  2004-6  Alberto Casagrande, Pieter Collins
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
 
/*! \file rational.h
 *  \brief Type definitions and conversion operators for rational numbers.
 */

#ifndef _ARIADNE_RATIONAL_H
#define _ARIADNE_RATIONAL_H

#include <cassert>

#include <gmpxx.h>

#include "../declarations.h"
#include "../numeric/numerical_traits.h"
#include "../numeric/arithmetic.h"
#include "../numeric/function.h"
#include "../numeric/integer.h"

namespace Ariadne {
  namespace Numeric {

#ifdef DOXYGEN

   /*!\ingroup Numeric
    * \brief A rational number.
    * 
    * An element of the field of rationals.
    * Must allow denotation of any rational.
    * May be created without loss of precision from any integral or floating point type, and from a dyadic.
    *
    * Currently implemented using mpq_class from the GNU Multiple Precision library.
    */
    class Rational { };
#else
    typedef mpq_class Rational;
#endif
  
    inline Integer numerator(const Rational& num){ 
      return num.get_num(); }
  
    inline Integer denominator(const Rational& num){ 
      return num.get_den();}
  
  
     
    template<> class numerical_traits<Rational> {
     public:
      typedef field_tag algebraic_category;
      typedef Rational field_extension_type;
    };
  
    template<> inline std::string name<Numeric::Rational>() { return "Rational"; }
    template<> inline std::string name<Numeric::Interval<Numeric::Rational> >() { return "Interval<Rational>"; }
    
    template<> inline double conv_approx(const Rational& x) { return x.get_d(); }
 
    template<> inline Rational conv_exact(const int& n) { return Rational(n); }
    template<> inline Rational conv_down(const int& n) { return conv_exact<Rational>(n); }
    template<> inline Rational conv_up(const int& n) { return conv_exact<Rational>(n); }
 
    template<> inline Rational conv_exact(const double& x) { return Rational(x); }
    template<> inline Rational conv_down(const double& x) { return conv_exact<Rational>(x); }
    template<> inline Rational conv_up(const double& x) { return conv_exact<Rational>(x); }
 
    template<> inline Rational conv_exact(const mpq_class& x) { return x; }
    template<> inline Rational conv_down(const mpq_class& x) { return conv_exact<Rational>(x); }
    template<> inline Rational conv_up(const mpq_class& x) { return conv_exact<Rational>(x); }
 
    template<> inline Rational neg_exact(const Rational& x) { return -x; }
    template<> inline Rational neg_down(const Rational& x) { return neg_exact(x); }
    template<> inline Rational neg_up(const Rational& x) { return neg_exact(x); }
    
    template<> inline Rational abs_exact(const Rational& x) { return (x>=0) ? x : Rational(-x); }
    template<> inline Rational abs_down(const Rational& x) { return abs_exact(x); }
    template<> inline Rational abs_up(const Rational& x) { return abs_exact(x); }
    
    template<> inline Rational add_exact(const Rational& x1, const Rational& x2) { return x1+x2; }
    template<> inline Rational add_down(const Rational& x1, const Rational& x2) { return add_exact(x1,x2); }
    template<> inline Rational add_up(const Rational& x1, const Rational& x2) { return add_exact(x1,x2); }
    template<> inline Rational add_approx(const Rational& x1, const Rational& x2) { return add_exact(x1,x2); }

    template<> inline Rational sub_exact(const Rational& x1, const Rational& x2) { return x1-x2; }
    template<> inline Rational sub_down(const Rational& x1, const Rational& x2) { return sub_exact(x1,x2); }
    template<> inline Rational sub_up(const Rational& x1, const Rational& x2) { return sub_exact(x1,x2); }
    template<> inline Rational sub_approx(const Rational& x1, const Rational& x2) { return add_exact(x1,x2); }
    
    template<> inline Rational mul_exact(const Rational& x1, const Rational& x2) { return x1*x2; }
    template<> inline Rational mul_down(const Rational& x1, const Rational& x2) { return mul_exact(x1,x2); }
    template<> inline Rational mul_up(const Rational& x1, const Rational& x2) { return mul_exact(x1,x2); }
    template<> inline Rational mul_approx(const Rational& x1, const Rational& x2) { return add_exact(x1,x2); }
    
    template<> inline Rational div_exact(const Rational& x1, const Rational& x2) { return x1/x2; }
    template<> inline Rational div_down(const Rational& x1, const Rational& x2) { return div_exact(x1,x2); }
    template<> inline Rational div_up(const Rational& x1, const Rational& x2) { return div_exact(x1,x2); }
    template<> inline  Rational div_approx(const Rational& x1, const Rational& x2) { return add_exact(x1,x2); }
  
    Rational quot(const Rational& x1, const Rational& x2);
  
  }

}

#endif /* _ARIADNE_RATIONAL_H */