/***************************************************************************
 *            function.decl.hpp
 *
 *  Copyright  2016-20  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/*! \file function.decl.hpp
 *  \brief Declarations of function types.
 */

#ifndef ARIADNE_FUNCTION_DECL_HPP
#define ARIADNE_FUNCTION_DECL_HPP

#include "../geometry/interval.decl.hpp"
#include "../geometry/box.decl.hpp"

namespace Ariadne {

template<class X> using Scalar=X;
template<class X> class Vector;

class Real;

//! \ingroup LinearAlgebraModule
//! \brief A scalar real number. Defined as a synomym for Real
using RealScalar=Scalar<Real>;
//! \ingroup LinearAlgebraModule
//! \brief A vector of real numbers. Defined as a synomym for Real
using RealVector=Vector<Real>;

// Domain declarations
typedef Interval<FloatDPValue> IntervalDomainType;
typedef Box<Interval<FloatDPValue>> BoxDomainType;

class RealDomain;
class EuclideanDomain;

// Function declarations
template<class P, class SIG> class Function;

//!@{
//! \name Type shorthands for function classes
template<class P, class... ARGS> using ScalarFunction = Function<P,RealScalar(ARGS...)>; //!< <p/> \ingroup FunctionModule
template<class P, class... ARGS> using VectorFunction = Function<P,RealVector(ARGS...)>; //!< <p/> \ingroup FunctionModule
template<class P, class RES> using UnivariateFunction = Function<P,RES(RealScalar)>; //!< <p/> \ingroup FunctionModule
template<class P, class RES> using MultivariateFunction = Function<P,RES(RealVector)>; //!< <p/> \ingroup FunctionModule

template<class P> using ScalarUnivariateFunction = Function<P,RealScalar(RealScalar)>; //!< <p/> \ingroup FunctionModule
template<class P> using VectorUnivariateFunction = Function<P,RealVector(RealScalar)>; //!< <p/> \ingroup FunctionModule
template<class P> using ScalarMultivariateFunction = Function<P,RealScalar(RealVector)>; //!< <p/> \ingroup FunctionModule
template<class P> using VectorMultivariateFunction = Function<P,RealVector(RealVector)>; //!< <p/> \ingroup FunctionModule

template<class SIG> using ApproximateFunction = Function<ApproximateTag,SIG>; //!< <p/> \ingroup FunctionModule
template<class SIG> using ValidatedFunction = Function<ValidatedTag,SIG>; //!< <p/> \ingroup FunctionModule
template<class SIG> using EffectiveFunction = Function<EffectiveTag,SIG>; //!< <p/> \ingroup FunctionModule
//!@}


//! \relates Function
//! \name Type synonyms for function classes
//!@{
using ApproximateScalarUnivariateFunction = ScalarUnivariateFunction<ApproximateTag>; //!< <p/>
using ValidatedScalarUnivariateFunction = ScalarUnivariateFunction<ValidatedTag>; //!< <p/>
using EffectiveScalarUnivariateFunction = ScalarUnivariateFunction<EffectiveTag>; //!< <p/>

using ApproximateScalarMultivariateFunction = ScalarMultivariateFunction<ApproximateTag>; //!< <p/>
using ValidatedScalarMultivariateFunction = ScalarMultivariateFunction<ValidatedTag>; //!< <p/>
using EffectiveScalarMultivariateFunction = ScalarMultivariateFunction<EffectiveTag>; //!< <p/>

using ApproximateVectorUnivariateFunction = VectorUnivariateFunction<ApproximateTag>; //!< <p/>
using ValidatedVectorUnivariateFunction = VectorUnivariateFunction<ValidatedTag>; //!< <p/>
using EffectiveVectorUnivariateFunction = VectorUnivariateFunction<EffectiveTag>; //!< <p/>

using ApproximateVectorMultivariateFunction = VectorMultivariateFunction<ApproximateTag>; //!< <p/>
using ValidatedVectorMultivariateFunction = VectorMultivariateFunction<ValidatedTag>; //!< <p/>
using EffectiveVectorMultivariateFunction = VectorMultivariateFunction<EffectiveTag>; //!< <p/>

using RealScalarUnivariateFunction = EffectiveScalarUnivariateFunction; //!< DEPRECATED
using RealScalarMultivariateFunction = EffectiveScalarMultivariateFunction; //!< DEPRECATED
using RealVectorUnivariateFunction = EffectiveVectorUnivariateFunction; //!< DEPRECATED
using RealVectorMultivariateFunction = EffectiveVectorMultivariateFunction; //!< DEPRECATED
//!@}

//! \ingroup FunctionModule
//! \brief An interface for generic functions.
//! \see Function
template<class P, class SIG> class FunctionInterface;

template<class P, class SIG> using FunctionPatch=Function<P,SIG>;
using ValidatedScalarUnivariateFunctionPatch = FunctionPatch<ValidatedTag,RealScalar(RealScalar)>;
using ValidatedVectorUnivariateFunctionPatch = FunctionPatch<ValidatedTag,RealVector(RealScalar)>;
using ValidatedScalarMultivariateFunctionPatch = FunctionPatch<ValidatedTag,RealScalar(RealVector)>;
using ValidatedVectorMultivariateFunctionPatch = FunctionPatch<ValidatedTag,RealVector(RealVector)>;


// Function models declarations


template<class P, class SIG, class PR=Void, class PRE=PR> class FunctionModelInterface;

template<class P, class SIG, class PR=Void, class PRE=PR> class FunctionModel;

//! \relates FunctionModel
//! \name Type shorthands
//!@{
template<class P, class ARG, class PR=Void, class PRE=PR> using ScalarFunctionModel = FunctionModel<P,RealScalar(ARG),PR,PRE>; //!< <p/>
template<class P, class ARG, class PR=Void, class PRE=PR> using VectorFunctionModel = FunctionModel<P,RealVector(ARG),PR,PRE>; //!< <p/>
template<class P, class RES, class PR=Void, class PRE=PR> using UnivariateFunctionModel = FunctionModel<P,RES(RealScalar),PR,PRE>; //!< <p/>
template<class P, class RES, class PR=Void, class PRE=PR> using MultivariateFunctionModel = FunctionModel<P,RES(RealVector),PR,PRE>; //!< <p/>
template<class P, class PR=Void, class PRE=PR> using ScalarUnivariateFunctionModel = FunctionModel<P,RealScalar(RealScalar),PR,PRE>; //!< <p/>
template<class P, class PR=Void, class PRE=PR> using VectorUnivariateFunctionModel = FunctionModel<P,RealVector(RealScalar),PR,PRE>; //!< <p/>
template<class P, class PR=Void, class PRE=PR> using ScalarMultivariateFunctionModel = FunctionModel<P,RealScalar(RealVector),PR,PRE>; //!< <p/>
template<class P, class PR=Void, class PRE=PR> using VectorMultivariateFunctionModel = FunctionModel<P,RealVector(RealVector),PR,PRE>; //!< <p/>

template<class PR=Void, class PRE=PR> using ValidatedScalarUnivariateFunctionModel = ScalarUnivariateFunctionModel<ValidatedTag,PR,PRE>; //!< <p/>
template<class PR=Void, class PRE=PR> using ValidatedVectorUnivariateFunctionModel = VectorUnivariateFunctionModel<ValidatedTag,PR,PRE>; //!< <p/>
template<class PR=Void, class PRE=PR> using ValidatedScalarMultivariateFunctionModel = ScalarMultivariateFunctionModel<ValidatedTag,PR,PRE>; //!< <p/>
template<class PR=Void, class PRE=PR> using ValidatedVectorMultivariateFunctionModel = VectorMultivariateFunctionModel<ValidatedTag,PR,PRE>; //!< <p/>

template<class PR=Void, class PRE=PR> using ApproximateScalarUnivariateFunctionModel = ScalarUnivariateFunctionModel<ApproximateTag,PR,PRE>; //!< <p/>
template<class PR=Void, class PRE=PR> using ApproximateVectorUnivariateFunctionModel = VectorUnivariateFunctionModel<ApproximateTag,PR,PRE>; //!< <p/>
template<class PR=Void, class PRE=PR> using ApproximateScalarMultivariateFunctionModel = ScalarMultivariateFunctionModel<ApproximateTag,PR,PRE>; //!< <p/>
template<class PR=Void, class PRE=PR> using ApproximateVectorMultivariateFunctionModel = VectorMultivariateFunctionModel<ApproximateTag,PR,PRE>; //!< <p/>

template<class P> using ScalarUnivariateFunctionModelDP = ScalarUnivariateFunctionModel<P,DoublePrecision>; //!< <p/>
template<class P> using VectorUnivariateFunctionModelDP = VectorUnivariateFunctionModel<P,DoublePrecision>; //!< <p/>
template<class P> using ScalarMultivariateFunctionModelDP = ScalarMultivariateFunctionModel<P,DoublePrecision>; //!< <p/>
template<class P> using VectorMultivariateFunctionModelDP = VectorMultivariateFunctionModel<P,DoublePrecision>; //!< <p/>

template<class P> using ScalarUnivariateFunctionModelMP = ScalarUnivariateFunctionModel<P,MultiplePrecision>; //!< <p/>
template<class P> using VectorUnivariateFunctionModelMP = VectorUnivariateFunctionModel<P,MultiplePrecision>; //!< <p/>
template<class P> using ScalarMultivariateFunctionModelMP = ScalarMultivariateFunctionModel<P,MultiplePrecision>; //!< <p/>
template<class P> using VectorMultivariateFunctionModelMP = VectorMultivariateFunctionModel<P,MultiplePrecision>; //!< <p/>
//!@}

//! \relates FunctionModel
//! \name Type synonyms
//!@{
using ValidatedScalarUnivariateFunctionModelDP = ScalarUnivariateFunctionModel<ValidatedTag,DoublePrecision>; //!< <p/>
using ValidatedVectorUnivariateFunctionModelDP = VectorUnivariateFunctionModel<ValidatedTag,DoublePrecision>; //!< <p/>
using ValidatedScalarMultivariateFunctionModelDP = ScalarMultivariateFunctionModel<ValidatedTag,DoublePrecision>; //!< <p/>
using ValidatedVectorMultivariateFunctionModelDP = VectorMultivariateFunctionModel<ValidatedTag,DoublePrecision>; //!< <p/>

using ApproximateScalarUnivariateFunctionModelDP = ScalarUnivariateFunctionModel<ApproximateTag,DoublePrecision>; //!< <p/>
using ApproximateVectorUnivariateFunctionModelDP = VectorUnivariateFunctionModel<ApproximateTag,DoublePrecision>; //!< <p/>
using ApproximateScalarMultivariateFunctionModelDP = ScalarMultivariateFunctionModel<ApproximateTag,DoublePrecision>; //!< <p/>
using ApproximateVectorMultivariateFunctionModelDP = VectorMultivariateFunctionModel<ApproximateTag,DoublePrecision>; //!< <p/>

using ValidatedScalarUnivariateFunctionModelMP = ScalarUnivariateFunctionModel<ValidatedTag,MultiplePrecision>; //!< <p/>
using ValidatedVectorUnivariateFunctionModelMP = VectorUnivariateFunctionModel<ValidatedTag,MultiplePrecision>; //!< <p/>
using ValidatedScalarMultivariateFunctionModelMP = ScalarMultivariateFunctionModel<ValidatedTag,MultiplePrecision>; //!< <p/>
using ValidatedVectorMultivariateFunctionModelMP = VectorMultivariateFunctionModel<ValidatedTag,MultiplePrecision>; //!< <p/>

using ApproximateScalarUnivariateFunctionModelMP = ScalarUnivariateFunctionModel<ApproximateTag,MultiplePrecision>; //!< <p/>
using ApproximateVectorUnivariateFunctionModelMP = VectorUnivariateFunctionModel<ApproximateTag,MultiplePrecision>; //!< <p/>
using ApproximateScalarMultivariateFunctionModelMP = ScalarMultivariateFunctionModel<ApproximateTag,MultiplePrecision>; //!< <p/>
using ApproximateVectorMultivariateFunctionModelMP = VectorMultivariateFunctionModel<ApproximateTag,MultiplePrecision>; //!< <p/>
//!@}

template<class P, class PR, class PRE=PR> struct FunctionModelTraits;

template<class F> class UnknownError;

template<class PR, class PRE> struct FunctionModelTraits<ValidatedTag,PR,PRE> {
    typedef RawFloat<PR> F; typedef RawFloat<PRE> FE;
    typedef Value<F> ValueType; typedef Error<FE> ErrorType;
    typedef PositiveUpperBound<F> NormType; typedef Interval<UpperBound<F>> RangeType;
    typedef Bounds<F> NumericType; typedef ValidatedNumber GenericNumericType;
    typedef F RawFloatType;
};
template<class PR> struct FunctionModelTraits<ApproximateTag,PR> {
    typedef RawFloat<PR> F;
    typedef Approximation<F> ValueType; typedef UnknownError<F> ErrorType;
    typedef PositiveApproximation<F> NormType; typedef Interval<Approximation<F>> RangeType;
    typedef Approximation<F> NumericType; typedef ApproximateNumber GenericNumericType;
    typedef F RawFloatType;
};

template<class P, class PR, class PRE=PR> using CanonicalNumericType = typename FunctionModelTraits<P,PR,PRE>::NumericType;
template<class P, class PR> using CanonicalCoefficientType = typename FunctionModelTraits<P,PR>::CoefficientType;
template<class P, class PRE> using CanonicalErrorType = typename FunctionModelTraits<P,PRE,PRE>::ErrorType;

template<class P> using CanonicalNumeric64Type = typename FunctionModelTraits<P,DoublePrecision>::NumericType;
template<class P> using CanonicalCoefficient64Type = typename FunctionModelTraits<P,DoublePrecision>::CoefficientType;
template<class P> using CanonicalError64Type = typename FunctionModelTraits<P,DoublePrecision,DoublePrecision>::ErrorType;

template<class X> using PrecisionType = typename X::PrecisionType;
template<class X> using ErrorPrecisionType = typename X::ErrorPrecisionType;



template<class P, class PR, class PRE=PR> class FunctionModelFactoryInterface;

template<class P, class PR, class PRE=PR> class FunctionModelFactory;
typedef FunctionModelFactory<ValidatedTag,DoublePrecision> ValidatedFunctionModelDPFactory;
template<class FMF, class D> class FunctionModelCreator;

class UniIndex;
class MultiIndex;

template<class I, class X> class Monomial;
template<class I, class X> class Polynomial;

//! \relates Monomial
//! \name Type shorthands and synonyms for Monomial classes.
//!@{
template<class X> using UnivariateMonomial = Monomial<UniIndex,X>; //!< <p/>
template<class X> using MultivariateMonomial = Monomial<MultiIndex,X>; //!< <p/>

using FloatDPApproximationMultivariateMonomial = MultivariateMonomial<FloatDPApproximation>; //!< <p/>
using FloatDPBoundsMultivariateMonomial = MultivariateMonomial<FloatDPBounds>; //!< <p/>
using FloatMPApproximationMultivariateMonomial = MultivariateMonomial<FloatMPApproximation>; //!< <p/>
using FloatMPBoundsMultivariateMonomial = MultivariateMonomial<FloatMPBounds>; //!< <p/>
//!@}

//! \relates Polynomial
//! \name Type shorthands and synonyms for Polynomial classes
//!@{
template<class X> using UnivariatePolynomial = Polynomial<UniIndex,X>; //!< <p/>
template<class X> using MultivariatePolynomial = Polynomial<MultiIndex,X>; //!< <p/>

using FloatDPApproximationMultivariatePolynomial = MultivariatePolynomial<FloatDPApproximation>; //!< <p/>
using FloatDPBoundsMultivariatePolynomial = MultivariatePolynomial<FloatDPBounds>; //!< <p/>
using FloatMPApproximationMultivariatePolynomial = MultivariatePolynomial<FloatMPApproximation>; //!< <p/>
using FloatMPBoundsMultivariatePolynomial = MultivariatePolynomial<FloatMPBounds>; //!< <p/>

using FloatDPApproximationMultivariatePolynomialVector = Vector<MultivariatePolynomial<FloatDPApproximation>>; //!< <p/>
using FloatDPBoundsMultivariatePolynomialVector = Vector<MultivariatePolynomial<FloatDPBounds>>; //!< <p/>
using FloatMPApproximationMultivariatePolynomialVector = Vector<MultivariatePolynomial<FloatMPApproximation>>; //!< <p/>
using FloatMPBoundsMultivariatePolynomialVector = Vector<MultivariatePolynomial<FloatMPBounds>>; //!< <p/>
//!@}


} // namespace Ariadne

#endif
