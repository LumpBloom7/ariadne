/***************************************************************************
 *            number.cpp
 *
 *  Copyright 2013--17  Pieter Collins
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

/*! \file number.cpp
 *  \brief
 */



#include "utility/module.hpp"
#include "numeric/paradigm.hpp"

#include "number.hpp"
#include "logical.hpp"
#include "integer.hpp"
#include "dyadic.hpp"
#include "rational.hpp"
#include "real.hpp"
#include "float64.hpp"
#include "floatmp.hpp"

#include "number_interface.hpp"
#include "number_wrapper.hpp"

/************ Number *********************************************************/

namespace Ariadne {

using NumberHandle = Handle<NumberInterface>;

// Declare operations for Integer and Rational numbers
Real sqrt(Real const&);
Real exp(Real const&);
Real log(Real const&);
Real sin(Real const&);
Real cos(Real const&);
Real tan(Real const&);
Real atan(Real const&);

template<> struct DispatchingTraits<Integer> { typedef Aware<Integer> AwareOfTypes; };
template class NumberWrapper<Integer>;
template<> struct DispatchingTraits<Dyadic> { typedef Aware<Dyadic,Integer> AwareOfTypes; };
template class NumberWrapper<Dyadic>;
template<> struct DispatchingTraits<Rational> { typedef Aware<Rational,Dyadic,Integer> AwareOfTypes; };
template class NumberWrapper<Rational>;

template<> struct DispatchingTraits<Real> { typedef Aware<Real> AwareOfTypes; };
template class NumberWrapper<Real>;

ExactDouble::operator ExactNumber() const { return Dyadic(*this).operator ExactNumber(); }
//ExactDouble::operator ExactNumber() const { if(std::isfinite(this->get_d()) { return Dyadic(*this).operator ExactNumber(); } else { return ExactNumber(new NumberWrapper<ExactDouble>(*this)); } }
Integer::operator ExactNumber() const { return ExactNumber(new NumberWrapper<Integer>(*this)); }
Dyadic::operator ExactNumber() const { return ExactNumber(new NumberWrapper<Dyadic>(*this)); }
Rational::operator ExactNumber() const { return ExactNumber(new NumberWrapper<Rational>(*this)); }
Real::operator EffectiveNumber() const { return EffectiveNumber(new NumberWrapper<Real>(*this)); }

template class NumberWrapper<Float64Approximation>;
//template class NumberWrapper<Float64LowerBound>;
//template class NumberWrapper<Float64UpperBound>;
template class NumberWrapper<Float64Bounds>;
template class NumberWrapper<Float64Ball>;
template class NumberWrapper<Float64Value>;

template<> Float64Approximation::operator ApproximateNumber() const { return ApproximateNumber(new NumberWrapper<Float64Approximation>(*this)); }
//template<> Float64LowerBound::operator ValidatedLowerNumber() const { return ValidatedLowerNumber(new NumberWrapper<Float64LowerBound>(*this)); }
//template<> Float64UpperBound::operator ValidatedUpperNumber() const { return ValidatedUpperNumber(new NumberWrapper<Float64UpperBound>(*this)); }
template<> Float64Bounds::operator ValidatedNumber() const { return ValidatedNumber(new NumberWrapper<Float64Bounds>(*this)); }
template<> Float64Ball::operator ValidatedNumber() const { return ValidatedNumber(new NumberWrapper<Float64Ball>(*this)); }
template<> Float64Value::operator ExactNumber() const { return ExactNumber(new NumberWrapper<Float64Value>(*this)); }

template class NumberWrapper<FloatMPApproximation>;
//template class NumberWrapper<FloatMPLowerBound>;
//template class NumberWrapper<FloatMPUpperBound>;
template class NumberWrapper<FloatMPBounds>;
template class NumberWrapper<FloatMPBall>;
template class NumberWrapper<FloatMPValue>;

template<> FloatMPApproximation::operator ApproximateNumber() const { return ApproximateNumber(new NumberWrapper<FloatMPApproximation>(*this)); }
//template<> FloatMPLowerBound::operator ValidatedLowerNumber() const { return ValidatedLowerNumber(new NumberWrapper<FloatMPLowerBound>(*this)); }
//template<> FloatMPUpperBound::operator ValidatedUpperNumber() const { return ValidatedUpperNumber(new NumberWrapper<FloatMPUpperBound>(*this)); }
template<> FloatMPBounds::operator ValidatedNumber() const { return ValidatedNumber(new NumberWrapper<FloatMPBounds>(*this)); }
template<> FloatMPBall::operator ValidatedNumber() const { return ValidatedNumber(new NumberWrapper<FloatMPBall>(*this)); }
template<> FloatMPValue::operator ExactNumber() const { return ExactNumber(new NumberWrapper<FloatMPValue>(*this)); }

template<> String class_name<NumberHandle>() { return "NumberHandle"; }

} // namespace Ariadne