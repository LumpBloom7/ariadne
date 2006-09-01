/***************************************************************************
 *            numerical_types.h
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
 
/*! \file numerical_types.h
 *  \brief Type definitions and conversion operators for fundamental Ariadne types.
 */

#ifndef _ARIADNE_NUMERICAL_TYPES_H
#define _ARIADNE_NUMERICAL_TYPES_H

#include <gmpxx.h>

#include "../declarations.h"
#include "../numeric/integer.h"
#include "../numeric/dyadic.h"
#include "../numeric/numerical_traits.h"
#include "../utility/stlio.h"

namespace Ariadne {
  namespace Numeric {

#ifdef DOXYGEN
    /*! \brief A 64-bit fixed-precision floating point number.
     *  \ingroup Numeric
     *
     * Standard operations are not exact, but must support interval arithmetic.
     *
     * Currently implemented by the built-in type double.
     */
    class Float64;
#else
    typedef double Float64;
#endif
     

#ifdef DOXYGEN
    /*! \brief A dyadic rational (i.e. of form \f$m/2^n\f$).
    *   \ingroup Numeric
    * 
    * A element of the ring of dyadic rationals.
    * Must allow denotation of any dyadic rational.
    * May be created without loss of precision from any integral or floating point type,
    * or from any rational of the form m/2^n.
    * May be created without loss of precision from any integral or floating point type,
    * or from any rational of the form m/2^n.
    *
    * Currently implemented using a modification of the Synaps dyadic class.
    *
    * FIXME: mpf_class does not implement addition, subtraction and multiplication exactly.
    */
    class Dyadic { };
#else
    typedef Synaps::dyadic Dyadic;
#endif
    

#ifdef DOXYGEN
    /*! \brief A multiple-precision floating-point type.
     *  \ingroup Numeric
     * 
     * Currently implemented using mpf_class from the GNU Multiple Precision library.
     */
    class MPFloat { };
#else
    typedef mpf_class MPFloat;
#endif
    

#ifdef DOXYGEN
    /*! \brief A rational number.
    *   \ingroup Numeric
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
  
  
    inline Integer precision(const MPFloat& num) {
      return mpf_get_prec(num.get_mpf_t());
    }
  
    inline Integer exponent(const MPFloat& num) {
      long int res; 
      mpf_get_d_2exp(&res,num.get_mpf_t()); 
      return res-1; 
    }
  
    inline MPFloat mantissa(const MPFloat& num) {
      long int exp; 
      mpf_class res; 
      mpf_get_d_2exp(&exp,num.get_mpf_t()); 
      exp=exp-1;
      mpf_div_2exp(res.get_mpf_t(),num.get_mpf_t(),exp); 
      return res; 
    }
  
    inline Dyadic mantissa(const Dyadic& num) {
      return num.mantissa();
    }
  
    inline int exponent(const Dyadic& num) {
      return num.exponent();
    }
  
    inline int precision(const Dyadic& num) {
      return num.precision();
    }
  
    inline Integer numerator(const Dyadic& num) {
      return num.numerator();
    }
  
    inline Integer denominator(const Dyadic& num) {
      return num.denominator();
    }
    
    inline Integer numerator(const Rational& num){ 
      return num.get_num(); }
  
    inline Integer denominator(const Rational& num){ 
      return num.get_den();}
  
  
    /* numerical traits */
    template<> class numerical_traits<Integer> {
     public:
      typedef field_tag algebraic_category;
      typedef Rational field_extension_type;
    };
  
    template<> class numerical_traits<double> {
     public:
      typedef field_tag algebraic_category;
      typedef double field_extension_type;
    };
  
    template<> class numerical_traits<MPFloat> {
     public:
      typedef ring_tag algebraic_category;
      typedef Rational field_extension_type;
    };
      
    template<> class numerical_traits<Dyadic> {
     public:
      typedef ring_tag algebraic_category;
      typedef Rational field_extension_type;
    };
      
    template<> class numerical_traits<Rational> {
     public:
      typedef field_tag algebraic_category;
      typedef Rational field_extension_type;
    };
  
  }
  
  namespace Base {

    template<> inline std::string name<Numeric::Float64>() { return "Float64"; }
    template<> inline std::string name<Numeric::MPFloat>() { return "MPFloat"; }
    template<> inline std::string name<Numeric::Dyadic>() { return "Dyadic"; }
    template<> inline std::string name<Numeric::Rational>() { return "Rational"; }
    template<> inline std::string name<Numeric::Interval<Numeric::Float64> >() { return "Interval<Float64>"; }
    template<> inline std::string name<Numeric::Interval<Numeric::MPFloat> >() { return "Interval<MPFloat>"; }
    template<> inline std::string name<Numeric::Interval<Numeric::Dyadic> >() { return "Interval<Dyadic>"; }
    template<> inline std::string name<Numeric::Interval<Numeric::Rational> >() { return "Interval<Rational>"; }
    
    template<typename R> inline R convert_to(const Numeric::MPFloat& x) 
    { return R(x); }
    template<typename R> inline R convert_to(const Numeric::Dyadic& x) 
    { return R(x); }
    template<typename R> inline R convert_to(const Numeric::Rational& x) 
    { return R(x); }
    template<typename R> inline R convert_to(const double& x) 
    { return R(x); } 

  }

}

#endif /* _ARIADNE_NUMERICAL_TYPES */