/***************************************************************************
 *            interval.h
 *
 *  Copyright 2008-10  Pieter Collins
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

/*! \file interval.h
 *  \brief ExactInterval numeric class.
 */

#ifndef ARIADNE_INTERVAL_H
#define ARIADNE_INTERVAL_H

#include <iostream>
#include <cassert>

#include "utility/declarations.h"

#include "utility/tribool.h"
#include "numeric/logical.h"
#include "numeric/rounding.h"
#include "numeric/real.h"
#include "numeric/float.h"

namespace Ariadne {

class Tribool;

class ExactInterval;
class UpperInterval;
class ApproximateInterval;

// Allow vector arithmetic on UpperInterval
template<class X> struct IsScalar;
template<> struct IsScalar<UpperInterval> : True { };

typedef ExactInterval ExactFloatInterval;
typedef UpperInterval UpperFloatInterval;
typedef ApproximateInterval ApproximateFloatInterval;

typedef UpperInterval NumericInterval;

class UnitInterval;

struct SizeOne { operator SizeType() const { return 1u; } };

//! \ingroup NumericModule
//! \brief Intervals with floating-point endpoints supporting outwardly-rounded arithmetic.
//! \details
//! Note that <c>%ExactInterval(3.3)</c> yields the singleton interval \f$[3.2999999999999998224,3.2999999999999998224]\f$ (the constant is first interpreted by the C++ compiler to give a C++ \c double, whereas <c>%ExactInterval("3.3")</c> yields the interval \f$[3.2999999999999998224,3.3000000000000002665]\f$ enclosing \f$3.3\f$.
//!
//! Comparison tests on \c ExactInterval use the idea that an interval represents a single number with an unknown value.
//! Hence the result is of type \c Tribool, which can take values { \c True, \c False, \c Indeterminate }.
//! Hence a test \f$[l_1,u_1]\leq [l_2,u_2]\f$ returns \c True if \f$u_1\leq u_2\f$, since in this case \f$x_1\leq x_2\f$ whenever \f$x_1\in[l_1,u_2]\f$ and \f$x_2\in[l_2,u_2]\f$, \c False if \f$l_1>u_2\f$, since in this case we know \f$x_1>x_2\f$, and \c Indeterminate otherwise, since in this case we can find \f$x_1,x_2\f$ making the result either true or false.
//! In the case of equality, the comparison \f$[l_1,u_1]\f$==\f$[l_2,u_2]\f$ only returns \c True if both intervals are singletons, since otherwise we can find values making the result either true of false.
//!
//! To obtain the lower and upper bounds of an interval, use \c ivl.lower() and \c ivl.upper_raw().
//! To obtain the midpoint and radius, use \c ivl.midpoint() and \c ivl.radius().
//! Alternatives \c midpoint(ivl) and \c radius(ivl) are also provided.
//! Note that \c midpoint and \c radius return approximations to the true midpoint and radius of the interval. If \f$m\f$ and \f$r\f$ are the returned midpoint and radius of the interval \f$[l,u]\f$, the using exact arithmetic, we guarentee \f$m-r\leq l\f$ and \f$m+r\geq u\f$
//!
//! To test if an interval contains a point or another interval, use \c encloses(ExactInterval,Float64) or \c encloses(ExactInterval,ExactInterval).
//! The test \c refines(ExactInterval,ExactInterval) can also be used.
//! \sa Float64
//!
//! \par Python interface
//!
//! In the Python interface, %Ariadne intervals can be constructed from Python literals of the form \c {a:b} or (deprecated) \c [a,b] .
//! The former is preferred, as it cannot be confused with literals for other classes such as Vector and Array types.
//! Automatic conversion is used to convert ExactInterval literals of the form \c {a,b} to an ExactInterval in functions.
//!
//! Care must be taken when defining intervals using floating-point coefficients, since values are first converted to the nearest
//! representable value by the Python interpreter. <br><br>
//! \code
//!   ExactInterval({1.1:2.3}) # Create the interval [1.1000000000000001, 2.2999999999999998]
//!   ExactInterval({2.5:4.25}) # Create the interval [2.5, 4.25], which can be represented exactly
//!   ExactInterval([2.5,4.25]) # Alternative syntax for creating the interval [2.5, 4.25]
//! \endcode
class ExactInterval {
  public:
    typedef Exact Paradigm;
    typedef ExactInterval NumericType;
  public:
    //! \brief Default constructor yields the singleton zero interval \a [0,0].
    ExactInterval() : l(0.0), u(0.0) { }
    template<class N, EnableIf<IsIntegral<N>> =dummy>
        ExactInterval(N n) : l(n), u(n) { }
    template<class NL, class NU, EnableIf<IsIntegral<NL>> =dummy, EnableIf<IsIntegral<NU>> =dummy>
        ExactInterval(NL nl, NU nu) : l(nl), u(nu) { }
    //! \brief Copy constructor.
    ExactInterval(const ExactInterval& i) : l(i.l), u(i.u) { }
    //! \brief Convert from a floating-point number with an exact representation.
    explicit ExactInterval(const ExactFloat64& x) : l(x.raw()), u(x.raw()) { }
    //! \brief Convert from a dyadic number.
    explicit ExactInterval(const Dyadic& x);
    //! \brief Convert from a decimal number.
    explicit ExactInterval(const Decimal& x);
    //! \brief Convert from a raw float number.
    explicit ExactInterval(const Float64& x);

    //! \brief Convert to a floating-point approximation.
    explicit operator Float64 () const { return half_exact(add_near(l.raw(),u.raw())); }
    //! \brief Convert from a floating-point number with an exact representation.
    explicit operator ValidatedFloat64 () const { return ValidatedFloat64(this->lower(),this->upper()); }

    //! \brief Create from explicitly given lower and upper bounds. Yields the interval \a [lower,upper].
    ExactInterval(const Float64& lower, const Float64& upper) : l(lower), u(upper) { }
    //! \brief Convert from a floating-point number with an exact representation.
    ExactInterval(const ExactFloat64& lower, const ExactFloat64& upper) : l(lower.raw()), u(upper.raw()) { }
    //! \brief Create from explicitly given lower and upper bounds. Yields the interval \a [lower,upper].
    ExactInterval(const Real& lower, const Real& upper);

    explicit ExactInterval(const Integer& z);
    explicit ExactInterval(const Rational& q);
    ExactInterval(const Rational& lower, const Rational& upper);

    SizeOne dimension() const { return SizeOne(); }

    //! \brief The lower bound of the interval.
    const Float64& lower_value() const { return l; }
    const Float64& lower_raw() const { return l; }
    //! \brief The upper bound of the interval.
    const Float64& upper_value() const { return u; }
    const Float64& upper_raw() const { return u; }

    //! \brief The lower bound of the interval.
    const ExactFloat64& lower() const { return reinterpret_cast<ExactFloat64 const&>(l); }
    //! \brief The upper bound of the interval.
    const ExactFloat64& upper() const { return reinterpret_cast<ExactFloat64 const&>(u); }
    //! \brief The midpoint of the interval.
    const ValidatedFloat64 midpoint() const { return half(this->lower()+this->upper()); }
    //! \brief The radius of the interval.
    const ValidatedFloat64 radius() const { return half(this->upper()-this->lower()); }
    //! \brief An over-approximation to the width of the interval.
    const ErrorFloat64 width() const { return ErrorFloat64(sub_up(u,l)); }

    //! \brief An approximation to the midpoint of the interval.
    const ExactFloat64 centre() const { return ExactFloat64(half(add_near(l,u))); }
    //! \brief An over-approximation to the distance between the centre and the endpoints.
    const ErrorFloat64 error() const { Float64 c=half(add_near(l,u)); return ErrorFloat64(max(sub_up(u,c),sub_up(c,l))); }

    //! \brief Tests if the interval is empty.
    Bool empty() const { return l>u; }
    //! \brief Tests if the interval is a singleton.
    Bool singleton() const { return l==u; }

    //! \brief Sets the interval to a "canonical" empty interval \a [1,0].
    Void set_empty() { l=+std::numeric_limits< double >::infinity(); u=-std::numeric_limits< double >::infinity(); }
    Void set_lower(const ExactFloat64& lower) { l=lower.raw(); }
        // ARIADNE_ASSERT(lower<=this->u);
    Void set_upper(const ExactFloat64& upper) { u=upper.raw(); }
        // ARIADNE_ASSERT(this->l<=upper);
    Void set(const ExactFloat64& lower, const ExactFloat64& upper) { l=lower.raw(); u=upper.raw(); }
        // ARIADNE_ASSERT(lower<=upper);
    Void set_lower(const Float64& lower) { l=lower; }
        // ARIADNE_ASSERT(lower<=this->u);
    Void set_upper(const Float64& upper) { u=upper; }
        // ARIADNE_ASSERT(this->l<=upper);
    Void set(const Float64& lower, const Float64& upper) { l=lower; u=upper; }
        // ARIADNE_ASSERT(lower<=upper);
  public:
    //! \brief Extract a double-precision point approximation to the value represented by the interval.
    double get_d() const { return (this->l.get_d()+this->u.get_d())/2; }
    static Nat output_precision;
    static Void set_output_precision(Nat p) { output_precision=p; }
  private:
    Float64 l, u;
};

OutputStream& operator<<(OutputStream& os, const ExactInterval& ivl);

inline ValidatedFloat64 midpoint(ExactInterval i) {
    return i.midpoint();
}

inline ValidatedFloat64 radius(ExactInterval i) {
    return i.radius();
}

inline ErrorFloat64 width(ExactInterval i) {
    return i.width();
}

//! \related ExactInterval \brief Test if the intervals are equal (as sets).
inline Bool equal(ExactInterval i1, ExactInterval i2) {
    //std::cerr<<"equal(i1,i2) with i1="<<i1<<"; i2="<<i2<<std::endl;
    return i1.lower_raw()==i2.lower_raw() && i1.upper_raw()==i2.upper_raw();
}

//! \related ExactInterval \brief Test if the interval is empty.
inline Bool empty(ExactInterval i) {
    return i.lower_raw()>i.upper_raw();
}

//! \related ExactInterval \brief Test if the interval is bounded.
inline Bool bounded(ExactInterval i) {
    return i.lower_raw()!=-inf.raw() && i.upper_raw()!=+inf.raw();
}

//! \related ExactInterval \brief The intersection of two intervals.
inline ExactInterval intersection(ExactInterval i1, ExactInterval i2) {
    return ExactInterval(max(i1.lower_raw(),i2.lower_raw()),min(i1.upper_raw(),i2.upper_raw()));
}

//! \related ExactInterval \brief The hull of two intervals, equal to the smallest interval containing both as subsets.
inline ExactInterval hull(ExactInterval i1, ExactInterval i2) {
    assert(i1.lower_raw()<=i1.upper_raw() && i2.lower_raw()<=i2.upper_raw());
    return ExactInterval(min(i1.lower_raw(),i2.lower_raw()),max(i1.upper_raw(),i2.upper_raw()));
}

//! \related ExactInterval \brief The hull of an interval and a point, equal to the smallest interval containing both.
inline ExactInterval hull(ExactInterval i1, ExactFloat64 x2) {
    return ExactInterval(min(i1.lower_raw(),x2.raw()),max(i1.upper_raw(),x2.raw()));
}

//! \related ExactInterval \brief The hull of two points, equal to the smallest interval containing both.
inline ExactInterval hull(ExactFloat64 x1, ExactFloat64 x2) {
    return ExactInterval(min(x1.raw(),x2.raw()),max(x1.raw(),x2.raw()));
}

inline Bool operator==(ExactInterval i1, ExactInterval i2) {
    return equal(i1,i2);
}

inline Bool operator!=(ExactInterval i1, ExactInterval i2) {
    return !equal(i1,i2);
}


//! \related ExactInterval \brief Test if the interval \a I contains the number \a x.
inline Bool contains(ExactInterval i, ExactFloat64 x) { return i.lower_raw()<=x.raw() && x.raw()<=i.upper_raw(); }
inline Bool contains(ExactInterval i, ValidatedFloat64 x) { return i.lower_raw()<=x.lower_raw() && x.upper_raw()<=i.upper_raw(); }
inline Bool contains(ExactInterval i, ApproximateFloat64 x) { return i.lower_raw()<=x.raw() && x.raw()<=i.upper_raw(); }
inline Bool contains(ExactInterval i, RawFloat64 x) { return i.lower_raw()<=x && x<=i.upper_raw(); }

inline Bool element(ExactFloat64 x, ExactInterval i) { return i.lower_raw()<=x.raw() && x.raw()<=i.upper_raw(); }
inline Bool element(ValidatedFloat64 x, ExactInterval i) { return i.lower_raw()<=x.lower_raw() && x.upper_raw()<=i.upper_raw(); }
inline Bool element(ApproximateFloat64 x, ExactInterval i) { return i.lower_raw()<=x.raw() && x.raw()<=i.upper_raw(); }
inline Bool element(RawFloat64 x, ExactInterval i) { return i.lower_raw()<=x && x<=i.upper_raw(); }

//! \related ExactInterval \brief Test if the interval \a I1 is a subset of \a I2.
inline Bool subset(ExactInterval i1, ExactInterval i2) { return i1.lower_raw()>=i2.lower_raw() && i1.upper_raw()<=i2.upper_raw(); }
//! \related ExactInterval \brief Test if the interval \a I1 is a superset of \a I2.
inline Bool superset(ExactInterval i1, ExactInterval i2) { return i1.lower_raw()<=i2.lower_raw() && i1.upper_raw()>=i2.upper_raw(); }
//! \related ExactInterval \brief Test if the interval \a I1 is disjoint from \a I2. Returns \c false even if the two intervals only have an endpoint in common.
inline Bool disjoint(ExactInterval i1, ExactInterval i2) { return i1.lower_raw()>i2.upper_raw() || i1.upper_raw()<i2.lower_raw(); }
//! \related ExactInterval \brief Test if the interval \a I1 intersects \a I2. Returns \c true even if the two intervals only have an endpoint in common.
inline Bool intersect(ExactInterval i1, ExactInterval i2) { return i1.lower_raw()<=i2.upper_raw() && i1.upper_raw()>=i2.lower_raw(); }

//! \related ExactInterval \brief Test if the closed interval \a I1 is disjoint from the closed interval \a I2.
//! Returns \c false if the two intervals only have an endpoint in common.
inline Bool separated(ExactInterval i1, ExactInterval i2) { return i1.lower_raw()>i2.upper_raw() || i1.upper_raw()<i2.lower_raw(); }
//! \related ExactInterval \brief Test if the interval \a I1 overlaps \a I2.
//! Returns \c false if the two intervals only have an endpoint in common.
//! Returns \c true if one of the intervals is a singleton in the interior of the other.
inline Bool overlap(ExactInterval i1, ExactInterval i2) { return i1.lower_raw()<i2.upper_raw() && i1.upper_raw()>i2.lower_raw(); }
//! \related ExactInterval \brief Test if the (closed) interval \a I1 is a subset of the interior of \a I2.
inline Bool inside(ExactInterval i1, ExactInterval i2) { return i1.lower_raw()>i2.lower_raw() && i1.upper_raw()<i2.upper_raw(); }
//! \related ExactInterval \brief Test if the interior of the interval \a I1 is a superset of the (closed) interval \a I2.
inline Bool covers(ExactInterval i1, ExactInterval i2) { return i1.lower_raw()<i2.lower_raw() && i1.upper_raw()>i2.upper_raw(); }


class UpperInterval;
inline ValidatedFloat64 make_singleton(UpperInterval const& ivl);
inline UpperInterval make_interval(ValidatedFloat64 const& x);

extern const UpperInterval pi_ivl;

//! \brief An over-approximation to an interval set.
class UpperInterval {
  public:
    typedef Upper Paradigm;
    typedef UpperInterval NumericType;

    //! \brief Default constructor yields the singleton zero interval \a [0,0].
    explicit UpperInterval() : l(0.0), u(0.0) { }
    //! \brief Construct a singleton interval.
    template<class F, EnableIf<IsSame<F,Float64>> =dummy> explicit UpperInterval(F point) : l(point), u(point) { }
    template<class N, EnableIf<IsIntegral<N>> = dummy> explicit UpperInterval(N point) : l(point), u(point) { }
    template<class D, EnableIf<IsFloatingPoint<D>> =dummy> explicit UpperInterval(D point) : l(point), u(point) { }
    //! \brief Create from explicitly given lower and upper bounds. Yields the interval \a [lower,upper].
    explicit UpperInterval(Float64 lower, Float64 upper) : l(lower), u(upper) { }
    template<class N, EnableIf<IsIntegral<N>> = dummy> UpperInterval(N lower, N upper) : l(lower), u(upper) { }

    //! \brief Construct an over-approximating interval.
    UpperInterval(LowerFloat64 lower, UpperFloat64 upper) : UpperInterval(lower.raw(),upper.raw()) { }
    //! \brief Convert from an exact interval.
    UpperInterval(ExactInterval ivl) : UpperInterval(ivl.lower_raw(),ivl.upper_raw()) { }

    //! \brief Construct a singleton interval.
    explicit UpperInterval(ExactFloat64 point) : l(point.raw()), u(point.raw()) { }
    explicit UpperInterval(ValidatedFloat64 point) : l(point.lower_raw()), u(point.upper_raw()) { }
    explicit UpperInterval(Real point) : UpperInterval(ValidatedFloat64(point)) { }

    UpperInterval& operator=(ValidatedFloat64 point) { l=point.lower_raw(); u=point.upper_raw(); }

    //! \brief Set the lower bound of the interval.
    Void set_lower(LowerFloat64 lower) { l=lower.raw(); }
    //! \brief Set the upper bound of the interval.
    Void set_upper(UpperFloat64 upper) { u=upper.raw(); }

    //! \brief The lower bound of the interval.
    const Float64& lower_raw() const { return l; }
    //! \brief The upper bound of the interval.
    const Float64& upper_raw() const { return u; }

    //! \brief The lower bound of the interval.
    const LowerFloat64& lower() const { return reinterpret_cast<LowerFloat64 const&>(l); }
    //! \brief The upper bound of the interval.
    const UpperFloat64& upper() const { return reinterpret_cast<UpperFloat64 const&>(u); }
    //! \brief The midpoint of the interval.
    const ApproximateFloat64 midpoint() const { return half(this->lower()+this->upper()); }
    const ExactFloat64 centre() const { return make_exact(half(this->lower()+this->upper())); }
    //! \brief An over-approximation to the width of the interval.
    const PositiveUpperFloat64 width() const { return PositiveUpperFloat64((this->upper()-this->lower()).raw()); }
    //! \brief The radius of the interval.
    const PositiveUpperFloat64 radius() const { return half(this->width()); }

    explicit operator ValidatedFloat64() const { return ValidatedFloat64(l,u); }

    Tribool empty() const { return Tribool(l>u) || Tribool(indeterminate); }

    friend const ApproximateFloat64 midpoint(UpperInterval const& ivl) {
        return ivl.midpoint(); }
    friend Bool empty(UpperInterval const& ivl) {
        return (ivl.l>ivl.u); }
    friend Bool bounded(UpperInterval const& ivl) {
        return -inf.raw()<ivl.lower_raw() && ivl.upper_raw()<+inf.raw(); }
    friend Bool contains(UpperInterval const& ivl, ExactFloat64 const& x) {
        return ivl.lower_raw() <= x.raw() && x.raw() <= ivl.upper_raw(); }
    friend Tribool inside(UpperInterval const& ivl1, ExactInterval const& ivl2) {
        return (ivl1.lower_raw()>ivl2.lower_raw() && ivl1.upper_raw()<ivl2.upper_raw()) || Tribool(indeterminate); }
    friend Tribool subset(UpperInterval const& ivl1, ExactInterval const& ivl2) {
        return (ivl1.l>=ivl2.lower_raw() && ivl1.u<=ivl2.upper_raw()) || Tribool(indeterminate); }
    friend Bool equal(UpperInterval const& ivl1, UpperInterval const& ivl2) {
        return ivl1.l==ivl2.l && ivl1.u==ivl2.u; }
    friend Bool refines(UpperInterval const& ivl1, UpperInterval const& ivl2) {
        return ivl1.l>=ivl2.l && ivl1.u<=ivl2.u; }
    friend Bool models(UpperInterval const& ivl1, ExactInterval const& ivl2) {
        return ivl1.l>=ivl2.lower_raw() && ivl1.u<=ivl2.upper_raw(); }
    friend Tribool disjoint(UpperInterval const& ivl1, UpperInterval const& ivl2) {
        return (ivl1.u<ivl2.l || ivl2.u<ivl1.l) || Tribool(indeterminate); }
    friend UpperInterval hull(UpperInterval const& ivl1, UpperInterval const& ivl2) {
        return UpperInterval(min(ivl1.l,ivl2.l),max(ivl1.u,ivl2.u)); }
    friend UpperInterval intersection(UpperInterval const& ivl1, UpperInterval const& ivl2) {
        return UpperInterval(max(ivl1.l,ivl2.l),min(ivl1.u,ivl2.u)); }
    friend UpperInterval widen(UpperInterval x) {
        return UpperInterval(sub_down(x.lower_raw(),Float64::min()),add_up(x.upper_raw(),Float64::max())); }
    friend OutputStream& operator<<(OutputStream& os, UpperInterval const& ivl) {
        return os << ExactInterval(ivl.lower_raw(),ivl.upper_raw()); }


    friend inline ValidatedFloat64 make_singleton(UpperInterval const& ivl) { return ValidatedFloat64(ivl.lower(),ivl.upper()); }
    friend inline UpperInterval make_interval(ValidatedFloat64 const& x) { return UpperInterval(x.lower(),x.upper()); }

    //! \related UpperInterval \brief The interval of possible maximum values. Yields the interval between \c i1.upper_raw() and \c i2.upper_raw().
    friend inline UpperInterval max(UpperInterval i1, UpperInterval i2) {
        return make_interval(max(make_singleton(i1),make_singleton(i2))); }
    //! \related UpperInterval \brief The interval of possible minimum values. Yields the interval between \c i1.lower() and \c i2.lower().
    friend inline UpperInterval min(UpperInterval i1, UpperInterval i2) {
        return make_interval(min(make_singleton(i1),make_singleton(i2))); }
    //! \related UpperInterval \brief The interval of possible absolute values. Yields \f$\{ |x| \mid x\in I\}\f$.
    friend inline UpperInterval abs(UpperInterval i) {
        return make_interval(abs(make_singleton(i))); }

    //! \related UpperInterval \brief Unary plus function. Yields the identity \f$I=\{+x | x\in I\}\f$.
    friend inline UpperInterval pos(UpperInterval i) {
        return make_interval(pos(make_singleton(i))); }
    //! \related UpperInterval \brief Unary negation function. Yields the exact interval \f$\{-x | x\in I\}\f$.
    friend inline UpperInterval neg(UpperInterval i) {
        return make_interval(neg(make_singleton(i))); }
    //! \related UpperInterval \brief Unary square function. Yields an over-approximation to \f$\{ x^2 \mid x\in I\}\f$.
    //! Note that if \a I contains positive and negative values, \c sqr(I) is tighter than \c I*I .
    friend UpperInterval sqr(UpperInterval i) {
        return make_interval(sqr(make_singleton(i))); }
    //! \related UpperInterval \brief Unary reciprocal function. Yields an over-approximation to \f$\{ 1/x \mid x\in I\}\f$.
    //! Yields \f$[-\infty,+\infty]\f$ if \a I contains \a 0 in its interior, and an interval containing \f$[1/u,+\infty]\f$ if \a I=[0,u] .
    friend UpperInterval rec(UpperInterval i) {
        return make_interval(rec(make_singleton(i))); }

    //! \related UpperInterval \brief Binary addition function. Yields an over-approximation to \f$\{ x_1+x_2 \mid x_1\in I_1 \wedge x_2\in I_2\}\f$.
    friend inline UpperInterval add(UpperInterval i1, UpperInterval i2) {
        return make_interval(add(make_singleton(i1),make_singleton(i2))); }
    //! \related UpperInterval \brief Subtraction function. Yields an over-approximation to \f$\{ x_1-x_2 \mid x_1\in I_1 \wedge x_2\in I_2\}\f$.
    friend inline UpperInterval sub(UpperInterval i1, UpperInterval i2) {
        return make_interval(sub(make_singleton(i1),make_singleton(i2))); }
    //! \related UpperInterval \brief Binary multiplication function. Yields an over-approximation to \f$\{ x_1\times x_2 \mid x_1\in I_1 \wedge x_2\in I_2\}\f$.
    friend UpperInterval mul(UpperInterval i1, UpperInterval i2) {
        return make_interval(mul(make_singleton(i1),make_singleton(i2))); }
    //! \related UpperInterval \brief Division function. Yields an over-approximation to \f$\{ x_1 \div x_2 \mid x_1\in I_1 \wedge x_2\in I_2\}\f$.
    friend UpperInterval div(UpperInterval i1, UpperInterval i2) {
        return make_interval(div(make_singleton(i1),make_singleton(i2))); }


    //! \related UpperInterval \brief Positive integer power function. Yields an over-approximation to \f$\{ x^m \mid x\in I\}\f$.
    friend UpperInterval pow(UpperInterval i, Nat m) {
        return make_interval(pow(make_singleton(i),m)); }
    //! \related UpperInterval \brief %Integer power function. Yields an over-approximation to \f$\{ x^n \mid x\in I\}\f$.
    friend UpperInterval pow(UpperInterval i, Int n) {
        return make_interval(pow(make_singleton(i),n)); }

    //! \related UpperInterval \brief Square-root function. Yields an over-approximation to \f$\{ \sqrt{x} \mid x\in I\}\f$.
    //! Requires \c I.lower()>=0 .
    friend UpperInterval sqrt(UpperInterval i) {
        return make_interval(sqrt(make_singleton(i))); }
    //! \related UpperInterval \brief Exponential function. Yields an over-approximation to \f$\{ \exp{x} \mid x\in I\}\f$.
    friend UpperInterval exp(UpperInterval i) {
        return make_interval(exp(make_singleton(i))); }
    //! \related UpperInterval \brief Natural logarithm function. Yields an over-approximation to \f$\{ \log{x} \mid x\in I\}\f$.
    //! Requires \c I.lower()>0 .
    friend UpperInterval log(UpperInterval i) {
        return make_interval(log(make_singleton(i))); }

    //! \related UpperInterval \brief Sine function. Yields an over-approximation to \f$\{ \sin{x} \mid x\in I\}\f$.
    friend UpperInterval sin(UpperInterval i) {
        return make_interval(sin(make_singleton(i))); }
    //! \related UpperInterval \brief Cosine function. Yields an over-approximation to \f$\{ \cos{x} \mid x\in I\}\f$.
    friend UpperInterval cos(UpperInterval i) {
        return make_interval(cos(make_singleton(i))); }
    //! \related UpperInterval \brief Tangent function. Yields an over-approximation to \f$\{ \tan{x} \mid x\in I\}\f$.
    friend UpperInterval tan(UpperInterval i) {
        return make_interval(tan(make_singleton(i))); }
    //! \related UpperInterval \brief Interse tangent function. Yields an over-approximation to \f$\{ \atan{x} \mid x\in I\}\f$.
    friend UpperInterval atan(UpperInterval i) {
        return make_interval(atan(make_singleton(i))); }

    //! \related UpperInterval \brief The magnitude of the interval \a I. Yields \f$ \max\{ |x|\,\mid\,x\in I \}\f$.
    friend inline PositiveUpperFloat64 mag(UpperInterval i) { return PositiveUpperFloat64(max(abs(i.lower_raw()),abs(i.upper_raw()))); }
    //! \related UpperInterval \brief The mignitude of the interval \a I. Yields \f$ \min\{ |x|\,\mid\,x\in I \}\f$.
    friend inline LowerFloat64 mig(UpperInterval i) { return LowerFloat64(min(Float64(0),min(abs(i.lower_raw()),abs(i.upper_raw())))); }

    //! \related UpperInterval \brief Unary plus operator. Should be implemented exactly and yield \f$\{ +x \mid x\in I\}\f$.
    friend inline UpperInterval operator+(const UpperInterval& i) { return pos(i); }
    //! \related UpperInterval \brief Unary negation operator. Should be implemented exactly and yield \f$\{ -x \mid x\in I\}\f$.
    friend inline UpperInterval operator-(const UpperInterval& i) { return neg(i); }
    //! \related UpperInterval \brief Binary addition operator. Guaranteed to yield an over approximation to \f$\{ x_1+x_2 \mid x_1\in I_1 \wedge x_2\in I_2\}\f$.
    friend inline UpperInterval operator+(const UpperInterval& i1, const UpperInterval& i2) { return add(i1,i2); }
    //! \related UpperInterval \brief Binary addition operator. Guaranteed to yield an over approximation to \f$\{ x_1-x_2 \mid x_1\in I_1 \wedge x_2\in I_2\}\f$.
    friend inline UpperInterval operator-(const UpperInterval& i1, const UpperInterval& i2) { return sub(i1,i2); }
    //! \related UpperInterval \brief Binary addition operator. Guaranteed to yield an over approximation to \f$\{ x_1*x_2 \mid x_1\in I_1 \wedge x_2\in I_2\}\f$.
    friend inline UpperInterval operator*(const UpperInterval& i1, const UpperInterval& i2) { return mul(i1,i2); }
    //! \related UpperInterval \brief Binary addition operator. Guaranteed to yield an over approximation to \f$\{ x_1/x_2 \mid x_1\in I_1 \wedge x_2\in I_2\}\f$. Yields \f$[-\infty,+\infty]\f$ if \f$0\in I_2\f$.
    friend inline UpperInterval operator/(const UpperInterval& i1, const UpperInterval& i2) { return div(i1,i2); };

    //! \related UpperInterval \brief Inplace addition operator.
    friend inline UpperInterval& operator+=(UpperInterval& i1, const UpperInterval& i2) { i1=add(i1,i2); return i1; }
    //! \related UpperInterval \brief Inplace subtraction operator.
    friend inline UpperInterval& operator-=(UpperInterval& i1, const UpperInterval& i2) { i1=sub(i1,i2); return i1; }
    //! \related UpperInterval \brief Inplace multiplication operator.
    friend inline UpperInterval& operator*=(UpperInterval& i1, const UpperInterval& i2) { i1=mul(i1,i2); return i1; }
    //! \related UpperInterval \brief Inplace division operator.
    friend inline UpperInterval& operator/=(UpperInterval& i1, const UpperInterval& i2) { i1=div(i1,i2); return i1; }

    // Standard equality operators
    //! \related UpperInterval \brief Equality operator. Tests equality of intervals as geometric objects, so \c [0,1]==[0,1] returns \c true.
    friend inline Bool operator==(const UpperInterval& i1, const UpperInterval& i2) { return i1.lower_raw()==i2.lower_raw() && i1.upper_raw()==i2.upper_raw(); }
    //! \related UpperInterval \brief Inequality operator. Tests equality of intervals as geometric objects, so \c [0,2]!=[1,3] returns \c true,
    //! even though the intervals possibly represent the same exact real value.
    friend inline Bool operator!=(const UpperInterval& i1, const UpperInterval& i2) { return i1.lower_raw()!=i2.lower_raw() || i1.upper_raw()!=i2.upper_raw(); }

    // Boost-style Tribool (in)equality operators
    //inline Tribool operator==(const UpperInterval& i1, const UpperInterval& i2) {
    //  if(i1.lower_raw()>i2.upper_raw() || i1.upper_raw()<i2.lower_raw()) { return false; } else if(i1.lower_raw()==i2.upper_raw() && i1.upper_raw()==i2.lower_raw()) { return true; } else { return indeterminate; } }
    //inline Tribool operator!=(const UpperInterval& i1, const UpperInterval& i2) { return !(i1==i2); }

    //! \related UpperInterval \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
    //! Hence \c [1.0,3.0]>[0.0,2.0] yields \c indeterminate since the first interval could represent the number 1.25 and the second 1.75.
    friend inline Tribool operator> (UpperInterval i1, UpperInterval i2) {
        if(i1.lower_raw()> i2.upper_raw()) { return true; }
        else if(i1.upper_raw()<=i2.lower_raw()) { return false; }
        else { return indeterminate; }
    }

    //! \related UpperInterval \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
    friend inline Tribool operator< (UpperInterval i1, UpperInterval i2) {
        if(i1.upper_raw()< i2.lower_raw()) { return true; }
        else if(i1.lower_raw()>=i2.upper_raw()) { return false; }
        else { return indeterminate; }
    }

    //! \related UpperInterval \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
    friend inline Tribool operator>=(UpperInterval i1, UpperInterval i2) {
        if(i1.lower_raw()>=i2.upper_raw()) { return true; }
        else if(i1.upper_raw()< i2.lower_raw()) { return false; }
        else { return indeterminate; }
    }

    //! \related UpperInterval \brief Strict greater-than comparison operator. Tests equality of represented real-point value.
    friend inline Tribool operator<=(UpperInterval i1, UpperInterval i2) {
        if(i1.upper_raw()<=i2.lower_raw()) { return true; }
        else if(i1.lower_raw()> i2.upper_raw()) { return false; }
        else { return indeterminate; }
    }

    // Mixed operations
    friend inline UpperInterval operator+(UpperInterval i1, ValidatedFloat64 x2) { return i1+make_interval(x2); }
    friend inline UpperInterval operator-(UpperInterval i1, ValidatedFloat64 x2) { return i1-make_interval(x2); }
    friend inline UpperInterval operator*(UpperInterval i1, ValidatedFloat64 x2) { return i1*make_interval(x2); }
    friend inline UpperInterval operator/(UpperInterval i1, ValidatedFloat64 x2) { return i1/make_interval(x2); }
    friend inline UpperInterval operator+(ValidatedFloat64 x1, UpperInterval i2) { return make_interval(x1)+i2; }
    friend inline UpperInterval operator-(ValidatedFloat64 x1, UpperInterval i2) { return make_interval(x1)-i2; }
    friend inline UpperInterval operator*(ValidatedFloat64 x1, UpperInterval i2) { return make_interval(x1)*i2; }
    friend inline UpperInterval operator/(ValidatedFloat64 x1, UpperInterval i2) { return make_interval(x1)/i2; }
    friend inline UpperInterval& operator+=(UpperInterval& i1, ValidatedFloat64 x2) { return i1+=make_interval(x2); }
    friend inline UpperInterval& operator-=(UpperInterval& i1, ValidatedFloat64 x2) { return i1-=make_interval(x2); }
    friend inline UpperInterval& operator*=(UpperInterval& i1, ValidatedFloat64 x2) { return i1*=make_interval(x2); }
    friend inline UpperInterval& operator/=(UpperInterval& i1, ValidatedFloat64 x2) { return i1/=make_interval(x2); }
    friend inline Tribool operator==(UpperInterval i1, ValidatedFloat64 x2) { return i1==make_interval(x2); }
    friend inline Tribool operator!=(UpperInterval i1, ValidatedFloat64 x2) { return i1!=make_interval(x2); }
    friend inline Tribool operator<=(UpperInterval i1, ValidatedFloat64 x2) { return i1<=make_interval(x2); }
    friend inline Tribool operator>=(UpperInterval i1, ValidatedFloat64 x2) { return i1>=make_interval(x2); }
    friend inline Tribool operator< (UpperInterval i1, ValidatedFloat64 x2) { return i1< make_interval(x2); }
    friend inline Tribool operator> (UpperInterval i1, ValidatedFloat64 x2) { return i1> make_interval(x2); }

  private:
    Float64 l, u;
};

template<class R, class A> inline R numeric_cast(A const&);
template<> inline UpperInterval numeric_cast(const Float64& a) { return UpperInterval(a,a); }
template<> inline Float64 numeric_cast(const UpperInterval& a) { return midpoint(a).raw(); }

const ApproximateFloat64 midpoint(UpperInterval const& ivl);
Bool bounded(UpperInterval const& ivl);
Bool contains(UpperInterval const& ivl, ExactNumber const& x);
Bool refines(UpperInterval const& ivl1, UpperInterval const& ivl2);
Bool models(UpperInterval const& ivl1, ExactInterval const& ivl2);
Tribool inside(UpperInterval const& ivl1, ExactInterval const& ivl2);
Tribool subset(UpperInterval const& ivl1, ExactInterval const& ivl2);
Tribool disjoint(UpperInterval const& ivl1, UpperInterval const& ivl2);
UpperInterval widen(UpperInterval i);

inline ExactInterval make_exact_interval(UpperInterval ivl) { return ExactInterval(ivl.lower_raw(),ivl.upper_raw()); }

// An interval one ulp wider
//! \related ExactInterval \brief An interval containing the given interval in its interior.
ExactInterval widen(ExactInterval i);
//! \related ExactInterval \brief An interval contained in the interior of the given interval.
ExactInterval narrow(ExactInterval i);

// Over-approximate by an interval with float coefficients
//! \related ExactInterval \brief Over-approximate the interval by one using builtin single-precision floating-point values as endpoints.
ExactInterval trunc(ExactInterval);
ExactInterval trunc(ExactInterval, Nat eps);

//! \related ExactInterval \brief The nearest representable number to the midpoint of the interval.
inline Float64 med(ExactInterval i) { return half_exact(add_near(i.lower_raw(),i.upper_raw())); }
//! \related ExactInterval \brief An over-approximation to the radius of the interval.
inline Float64 rad(ExactInterval i) { return half_exact(sub_up(i.upper_raw(),i.lower_raw())); }
//! \related ExactInterval \brief An over-approximation to the width of the interval.
inline ExactFloat64 diam(ExactInterval i) { return ExactFloat64(sub_up(i.upper_raw(),i.lower_raw())); }

inline UpperInterval hull(UpperInterval i1, ValidatedFloat64 x2) { return hull(i1,make_interval(x2)); }
inline UpperInterval hull(ValidatedFloat64 x1, UpperInterval i2) { return hull(make_interval(x1),i2); }
inline UpperInterval hull(ValidatedFloat64 x1, ValidatedFloat64 x2) { return hull(make_interval(x1),make_interval(x2)); }


#ifdef ARIADNE_ENABLE_SERIALIZATION
  template<class A> Void serialize(A& a, ExactInterval& ivl, const Nat version) {
    a & ivl.lower_raw() & ivl.upper_raw(); }
#endif

OutputStream& operator<<(OutputStream&, const ExactInterval&);
InputStream& operator>>(InputStream&, ExactInterval&);


//! \brief An over-approximation to an interval set.
class ApproximateInterval {
  public:
    typedef Approximate Paradigm;
    explicit ApproximateInterval() : l(0.0), u(0.0) { }
    explicit ApproximateInterval(Float64 point) : l(point), u(point) { }
    explicit ApproximateInterval(Float64 lower, Float64 upper) : l(lower), u(upper) { }
    explicit ApproximateInterval(ApproximateFloat64 point) : l(point.raw()), u(point.raw()) { }
    explicit ApproximateInterval(ApproximateFloat64 lower, ApproximateFloat64 upper) : l(lower.raw()), u(upper.raw()) { }
    ApproximateInterval(ExactInterval ivl) : l(ivl.lower_raw()), u(ivl.upper_raw()) { }
    ApproximateInterval(UpperInterval ivl) : l(ivl.lower_raw()), u(ivl.upper_raw()) { }
    Float64 const& lower_raw() const { return l; }
    Float64 const& upper_raw() const { return l; }
    ApproximateFloat64 lower() const { return ApproximateFloat64(l); }
    ApproximateFloat64 upper() const { return ApproximateFloat64(u); }
    ApproximateFloat64 midpoint() const { return ApproximateFloat64((l+u)/2); }
    ApproximateFloat64 radius() const { return ApproximateFloat64((u-l)/2); }
    ApproximateFloat64 width() const { return ApproximateFloat64(u-l); }
    friend Bool contains(ApproximateInterval const& ivl, ApproximateFloat64 const& x) {
        return ivl.lower_raw()<=x.raw() && x.raw()<=ivl.upper_raw(); }
    friend OutputStream& operator<<(OutputStream& os, const ApproximateInterval& ivl) {
        return os << ExactInterval(ivl.lower_raw(),ivl.upper_raw()); }
  private:
    Float64 l, u;
};


class UnitInterval
    : public ExactInterval
{
  public:
    UnitInterval() : ExactInterval(-1,+1) { }
};

} // namespace Ariadne

#endif
