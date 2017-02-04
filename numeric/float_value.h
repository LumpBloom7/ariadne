/***************************************************************************
 *            float_value.h
 *
 *  Copyright 2008-16  Pieter Collins
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

/*! \file float_value.h
 *  \brief Exact Floating-point representations of real numbers.
 */

#ifndef ARIADNE_FLOAT_VALUE_H
#define ARIADNE_FLOAT_VALUE_H

#include "utility/macros.h"

#include "number.decl.h"
#include "float.decl.h"

#include "logical.h"
#include "builtin.h"
#include "twoexp.h"

namespace Ariadne {

template<class PR> struct NumericTraits<FloatValue<PR>> {
    typedef ExactNumber GenericType;
    typedef PositiveFloatValue<PR> PositiveType;
    typedef Boolean LessType;
    typedef Boolean EqualsType;
};

static_assert(not IsGenericNumericType<FloatValue<Precision64>>::value,"");
static_assert(not IsGenericNumericType<FloatValue<PrecisionMP>>::value,"");

//! \ingroup NumericModule
//! \related Float64, Float64Bounds
//! \brief A floating-point number, which is taken to represent the \em exact value of a real quantity.
template<class PR> class FloatValue
    : DispatchNumericOperations<FloatValue<PR>,FloatBounds<PR>>
    , DispatchComparisonOperations<FloatValue<PR>,Boolean>
    , DefineMixedComparisonOperators<FloatValue<PR>,ExactNumber,Boolean>
    , DefineMixedComparisonOperators<FloatValue<PR>,Rational,Boolean>
//    , DefineMixedComparisonOperators<FloatValue<PR>,Dyadic,Boolean>
//    , DefineMixedComparisonOperators<FloatValue<PR>,Integer,Boolean>
//    , DefineMixedComparisonOperators<FloatValue<PR>,Int,Boolean>
//        , public DispatchFloatOperations<FloatBall<PR>>
        , public DispatchFloatOperations<FloatBounds<PR>>
    , DefineConcreteGenericArithmeticOperators<FloatValue<PR>>
    , DefineConcreteGenericComparisonOperators<FloatValue<PR>>
{
    typedef ExactTag P; typedef RawFloat<PR> FLT;
  public:
    typedef ExactTag Paradigm;
    typedef FloatValue<PR> NumericType;
    typedef ExactNumber GenericType;
    typedef FLT RawFloatType;
    typedef PR PrecisionType;
    typedef PR PropertiesType;
  public:
    FloatValue<PR>() : _v(0.0) { }
    explicit FloatValue<PR>(PrecisionType pr) : _v(0.0,pr) { }
    explicit FloatValue<PR>(RawFloatType const& v) : _v(v) { }

    template<class N, EnableIf<IsBuiltinIntegral<N>> = dummy> FloatValue<PR>(N n, PR pr) : FloatValue<PR>(ExactDouble(n),pr) { }
    FloatValue<PR>(ExactDouble d, PR pr);
    FloatValue<PR>(const Integer& z, PR pr);
    FloatValue<PR>(const TwoExp& t, PR pr);
    FloatValue<PR>(const Dyadic& w, PR pr);
    FloatValue<PR>(const FloatValue<PR>& x, PR pr);

    template<class N, EnableIf<IsBuiltinIntegral<N>> = dummy> FloatValue<PR>& operator=(N n) { _v=n; return *this; }
    FloatValue<PR>& operator=(const Integer& z);
    FloatValue<PR>& operator=(const TwoExp& t);
    FloatValue<PR>& operator=(const Dyadic& w);

    operator ExactNumber () const;
    explicit operator Dyadic () const;
    explicit operator Rational () const;

    FloatBall<PR> create(ValidatedNumber const&) const;
//    explicit operator RawFloatType () const { return _v; }

    PrecisionType precision() const { return _v.precision(); }
    PropertiesType properties() const { return _v.precision(); }
    GenericType generic() const { return this->operator GenericType(); }
    RawFloatType const& raw() const { return _v; }
    RawFloatType& raw() { return _v; }
    double get_d() const { return _v.get_d(); }

    FloatBall<PR> pm(FloatError<PR> _e) const;
  public:
    friend FloatValue<PR> operator*(FloatValue<PR> const&, TwoExp const&);
    friend FloatValue<PR> operator/(FloatValue<PR> const&, TwoExp const&);
    friend FloatValue<PR>& operator*=(FloatValue<PR>&, TwoExp const&);
    friend FloatValue<PR>& operator/=(FloatValue<PR>&, TwoExp const&);
    friend FloatError<PR> mag(FloatValue<PR> const&);
    friend FloatLowerBound<PR> mig(FloatValue<PR> const&);
    friend Bool same(FloatValue<PR> const&, FloatValue<PR> const&);
    friend OutputStream& operator<<(OutputStream&, FloatValue<PR> const&);
  public:
    friend Comparison cmp(FloatValue<PR> const& x1, Rational const& q2) { return cmp(x1.raw(),q2); }
    friend Comparison cmp(FloatValue<PR> const& x1, Dyadic const& w2) { return cmp(x1.raw(),w2); }
    friend Comparison cmp(FloatValue<PR> const& x1, Integer const& z2) { return cmp(x1.raw(),z2); }
  public:
    static Nat output_places;
    static Void set_output_places(Nat p) { output_places=p; }
  private: public:
    RawFloatType _v;
  private:
    friend FloatValue<PR> operator*(FloatValue<PR> const& x, TwoExp const& y) {
        return FloatValue<PR>(x.raw()*RawFloat<PR>(y,x.precision())); }
    friend FloatValue<PR> operator/(FloatValue<PR> const& x, TwoExp const& y) {
        return FloatValue<PR>(x.raw()/RawFloat<PR>(y,x.precision())); }
    friend FloatValue<PR>& operator*=(FloatValue<PR>& x, TwoExp const& y) { return x=x*y; }
    friend FloatValue<PR>& operator/=(FloatValue<PR>& x, TwoExp const& y) { return x=x/y; }
    friend OutputStream& operator<<(OutputStream& os, FloatValue<PR> const& x) {
        return Operations<FloatValue<PR>>::_write(os,x); }
};

template<class PR> class Positive<FloatValue<PR>> : public FloatValue<PR> {
  public:
    Positive<FloatValue<PR>>() : FloatValue<PR>() { }
    explicit Positive<FloatValue<PR>>(PR const& pr) : FloatValue<PR>(pr) { }
    template<class M, EnableIf<IsBuiltinUnsignedIntegral<M>> =dummy> Positive<FloatValue<PR>>(M m, PR pr) : FloatValue<PR>(m,pr) { }
    Positive<FloatValue<PR>>(TwoExp const& ex, PR pr) : FloatValue<PR>(ex,pr) { }
    explicit Positive<FloatValue<PR>>(Dyadic const& w, PR pr) : FloatValue<PR>(w,pr) { }
    explicit Positive<FloatValue<PR>>(RawFloat<PR> const& x) : FloatValue<PR>(x) { }
    explicit Positive<FloatValue<PR>>(FloatValue<PR> const& x) : FloatValue<PR>(x) { }
  public:
    friend Positive<FloatValue<PR>> hlf(Positive<FloatValue<PR>> const&);
};

template<class PR> inline PositiveFloatValue<PR> cast_positive(FloatValue<PR> const& x) {
    return PositiveFloatValue<PR>(x); }

static_assert(IsSame<decltype(declval<Float64Value>() < declval<Rational>()),Boolean>::value,"");

}

#endif
