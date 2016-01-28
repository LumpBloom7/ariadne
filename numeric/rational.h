/***************************************************************************
 *            numeric/rational.h
 *
 *  Copyright 2013-14  Pieter Collins
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

/*! \file numeric/rational.h
 *  \brief
 */



#ifndef ARIADNE_RATIONAL_H
#define ARIADNE_RATIONAL_H

#include "external/gmp.h"
#include "utility/typedefs.h"
#include "utility/metaprogramming.h"
#include "utility/string.h"
#include "numeric/integer.h"
#include "numeric/arithmetic.h"
#include "numeric/number.decl.h"
#include "numeric/float.decl.h"

namespace Ariadne {

class Int64;
class Float64;

enum class Comparison : char;


/************ Rational *******************************************************/

//! \ingroup UserNumericTypeSubModule
//! \brief %Rational numbers.
class Rational
    : Field<Rational>
    , DirectedLattice<Rational>
    , Ordered<Rational,Boolean>
    , DefineArithmeticOperators<Rational>
    , DefineComparisonOperators<Rational,Boolean>
{
  public:
    mpq_t _mpq;
  public:
    typedef ExactTag Paradigm;
    typedef Rational NumericType;
  public:
    ~Rational();
    Rational();
    Rational(const Integer&, const Integer&);
    template<class N, EnableIf<IsIntegral<N>> = dummy> Rational(N n);
    Rational(Int64);
    explicit Rational(Float64 const&);
    Rational(const Integer&);
    Rational(const Dyadic&);
    explicit Rational(const String&);
    explicit Rational(const Float64Value&);
    explicit Rational(const mpq_t);
    Rational(const Rational&);
    Rational(Rational&&);
    Rational& operator=(const Rational&);
    Rational& operator=(Rational&&);
    operator Number<ExactTag> () const;
    Integer get_num() const;
    Integer get_den() const;
    Integer numerator() const;
    Natural denominator() const;
    friend Rational operator/(Integer const& z1, Integer const& z2);

    friend Comparison cmp(Rational const& q1, Rational const& q2);

    friend OutputStream& operator<<(OutputStream& os, Rational const& q);
    friend InputStream& operator>>(InputStream& os, Rational& q);
    friend Rational operator"" _q(long double x);
  public:
    double get_d() const;
    mpq_t const& get_mpq() const;
  private:
    friend class Dyadic;
  private:
    explicit Rational(double, std::nullptr_t dummy);
};
template<> struct IsNumericType<Rational> : True { };
Rational operator"" _q(long double x);

template<class N, EnableIf<IsIntegral<N>>> inline Rational::Rational(N n) : Rational(Int64(n)) { }



} // namespace Ariadne

#endif
