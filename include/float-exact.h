/***************************************************************************
 *            float-exact.h
 *
 *  Copyright 2008-14  Pieter Collins
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

/*! \file float-exact.h
 *  \brief Exact floating-point number class, a subset of dyadic numbers.
 */
#ifndef ARIADNE_FLOAT_EXACT_H
#define ARIADNE_FLOAT_EXACT_H

#include <iostream> // For std::floor std::ceil etc
#include <iomanip> // For std::setprecision
#include <cmath> // For std::floor std::ceil etc
#include <algorithm> // For std::max, std::min
#include <limits> // For std::numeric_limits<double>

#include "rounding.h"
#include "rational.h"
#include "float.h"

namespace Ariadne {

//! \ingroup NumericModule
//! \related Float, ValidatedFloat
//! \brief A floating-point number, which is taken to represent the \em exact value of a real quantity.
class ExactFloat {
  private:
    Float _v;
  public:
    //! \brief Default constructor creates the number 0 (zero).
    ExactFloat() : _v(0) { }
    //! \brief Convert from a built-in positive integer.
    ExactFloat(unsigned int n) : _v(n) { }
    //! \brief Convert from a built-in integer.
    ExactFloat(int n) : _v(n) { }
    //! \brief Explicit construction from a built-in double-precision value.
    //! \details Tests to ensure that the number is not 'accidentally' created from a rounded version of a string literal,
    //! by comparing the input with it's single-precision approximation.
    explicit ExactFloat(double x) : _v(x) { }
    //! \brief Explicit construction from an approximate floating-point value.
    explicit ExactFloat(const Float& x) : _v(x) { }
#ifdef HAVE_GMPXX_H
    //! \brief Convert to a rational number.
    explicit operator Rational () const;
#endif
    //! \brief Explicit conversion to raw data type.
    explicit operator Float () const { return _v; }
    //! \brief The raw floating-point number with the same value.
    Float const& value() const { return _v; }
    //! \brief A double-precision approximateion.
    double get_d() const { return _v.get_d(); }

  public:
    static uint output_precision;
    static void set_output_precision(uint p) { output_precision=p; }
};

class ValidatedFloat;

inline ExactFloat operator+(const ExactFloat& x) { return ExactFloat(pos_exact(x.value())); }
inline ExactFloat operator-(const ExactFloat& x) { return ExactFloat(neg_exact(x.value())); }
inline ValidatedFloat operator+(const ExactFloat& x1,  const ExactFloat& x2);
inline ValidatedFloat operator-(const ExactFloat& x1,  const ExactFloat& x2);
inline ValidatedFloat operator*(const ExactFloat& x1,  const ExactFloat& x2);
inline ValidatedFloat operator/(const ExactFloat& x1,  const ExactFloat& x2);
inline ValidatedFloat operator+(const ValidatedFloat& x1,  const ExactFloat& x2);
inline ValidatedFloat operator-(const ValidatedFloat& x1,  const ExactFloat& x2);
inline ValidatedFloat operator*(const ValidatedFloat& x1,  const ExactFloat& x2);
inline ValidatedFloat operator/(const ValidatedFloat& x1,  const ExactFloat& x2);
inline ValidatedFloat operator+(const ExactFloat& x1,  const ValidatedFloat& x2);
inline ValidatedFloat operator-(const ExactFloat& x1,  const ValidatedFloat& x2);
inline ValidatedFloat operator*(const ExactFloat& x1,  const ValidatedFloat& x2);
inline ValidatedFloat operator/(const ExactFloat& x1,  const ValidatedFloat& x2);
inline ValidatedFloat pow(const ExactFloat& x, int n);
inline ValidatedFloat operator/(int n1,  const ExactFloat& x2);

inline std::ostream& operator<<(std::ostream& os, const ExactFloat& x) {
    return os << std::showpoint << std::setprecision(ExactFloat::output_precision) << x.value(); }

inline bool operator==(const ExactFloat& x1, const ExactFloat& x2) { return x1.value()==x2.value(); }
inline bool operator!=(const ExactFloat& x1, const ExactFloat& x2) { return x1.value()!=x2.value(); }
inline bool operator<=(const ExactFloat& x1, const ExactFloat& x2) { return x1.value()<=x2.value(); }
inline bool operator>=(const ExactFloat& x1, const ExactFloat& x2) { return x1.value()>=x2.value(); }
inline bool operator< (const ExactFloat& x1, const ExactFloat& x2) { return x1.value()< x2.value(); }
inline bool operator> (const ExactFloat& x1, const ExactFloat& x2) { return x1.value()> x2.value(); }

inline bool operator==(const ExactFloat& x1, double x2) { return x1.value()==x2; }
inline bool operator!=(const ExactFloat& x1, double x2) { return x1.value()!=x2; }
inline bool operator<=(const ExactFloat& x1, double x2) { return x1.value()<=x2; }
inline bool operator>=(const ExactFloat& x1, double x2) { return x1.value()>=x2; }
inline bool operator< (const ExactFloat& x1, double x2) { return x1.value()< x2; }
inline bool operator> (const ExactFloat& x1, double x2) { return x1.value()> x2; }

class ApproximateFloat;
inline const ExactFloat& make_exact(const Float& x) { return reinterpret_cast<const ExactFloat&>(x); }
inline const ExactFloat& make_exact(const ApproximateFloat& x) { return reinterpret_cast<const ExactFloat&>(x); }
template<template<class>class T> inline const T<ExactFloat>& make_exact(const T<ApproximateFloat>& t) { return reinterpret_cast<const T<ExactFloat>&>(t); }
template<template<class>class T> inline const T<ExactFloat>& make_exact(const T<Float>& t) { return reinterpret_cast<const T<ExactFloat>&>(t); }

//! \related Float \brief The constant infinity
//extern ExactFloat inf;

inline ExactFloat half(ExactFloat x) { return ExactFloat(half_exact(x.value())); }
inline ExactFloat pos(ExactFloat x) { return ExactFloat(pos_exact(x.value())); }
inline ExactFloat neg(ExactFloat x) { return ExactFloat(neg_exact(x.value())); }

inline ValidatedFloat sqr(ExactFloat x);
inline ValidatedFloat rec(ExactFloat x);
inline ValidatedFloat add(ExactFloat x, ExactFloat y);
inline ValidatedFloat sub(ExactFloat x, ExactFloat y);
inline ValidatedFloat mul(ExactFloat x, ExactFloat y);
inline ValidatedFloat div(ExactFloat x, ExactFloat y);
inline ValidatedFloat pow(ExactFloat x, int n);

inline ValidatedFloat rad(ExactFloat x, ExactFloat y);
inline ValidatedFloat med(ExactFloat x, ExactFloat y);


#ifdef HAVE_GMPXX_H
inline bool operator==(const ExactFloat& x, const Rational& q) { return x.get_d()==static_cast<const mpq_class&>(q); }
inline bool operator!=(const ExactFloat& x, const Rational& q) { return x.get_d()!=static_cast<const mpq_class&>(q); }
inline bool operator<=(const ExactFloat& x, const Rational& q) { return x.get_d()<=static_cast<const mpq_class&>(q); }
inline bool operator>=(const ExactFloat& x, const Rational& q) { return x.get_d()>=static_cast<const mpq_class&>(q); }
inline bool operator< (const ExactFloat& x, const Rational& q) { return x.get_d()< static_cast<const mpq_class&>(q); }
inline bool operator> (const ExactFloat& x, const Rational& q) { return x.get_d()> static_cast<const mpq_class&>(q); }

inline bool operator==(const Rational& q, const ExactFloat& x) { return static_cast<mpq_class>(q)==x.get_d(); }
inline bool operator!=(const Rational& q, const ExactFloat& x) { return static_cast<mpq_class>(q)!=x.get_d(); }
inline bool operator<=(const Rational& q, const ExactFloat& x) { return static_cast<mpq_class>(q)<=x.get_d(); }
inline bool operator>=(const Rational& q, const ExactFloat& x) { return static_cast<mpq_class>(q)>=x.get_d(); }
inline bool operator< (const Rational& q, const ExactFloat& x) { return static_cast<mpq_class>(q)< x.get_d(); }
inline bool operator> (const Rational& q, const ExactFloat& x) { return static_cast<mpq_class>(q)> x.get_d(); }
#endif // HAVE_GMPXX_H




} // namespace Ariadne

#endif