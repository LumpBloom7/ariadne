/***************************************************************************
 *            function/function_mixin.hpp
 *
 *  Copyright  2008-20  Pieter Collins
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

#ifndef ARIADNE_FUNCTION_MIXIN_HPP
#define ARIADNE_FUNCTION_MIXIN_HPP

#include "function/function_interface.hpp"

// Adaptors for classes to conform to the Function interface.

namespace Ariadne {

typedef Differential<FloatDPApproximation> ApproximateDifferentialDP;
typedef Differential<FloatDPBounds> ValidatedDifferentialDP;
typedef UnivariateDifferential<FloatDPApproximation> ApproximateUnivariateDifferentialDP;
typedef UnivariateDifferential<FloatDPBounds> ValidatedUnivariateDifferentialDP;
typedef TaylorModel<ApproximateTag,FloatDP> ApproximateTaylorModelDP;
typedef TaylorModel<ValidatedTag,FloatDP> ValidatedTaylorModelDP;
typedef TaylorModel<ApproximateTag,FloatMP> ApproximateTaylorModelMP;
typedef TaylorModel<ValidatedTag,FloatMP> ValidatedTaylorModelMP;

typedef TaylorModel<ValidatedTag,FloatDPUpperInterval> ValidatedIntervalTaylorModelDP;

typedef Formula<ApproximateNumber> ApproximateFormula;
typedef Formula<ValidatedNumber> ValidatedFormula;
typedef Formula<EffectiveNumber> EffectiveFormula;
typedef Algebra<ApproximateNumber> ApproximateAlgebra;
typedef Algebra<ValidatedNumber> ValidatedAlgebra;
typedef Algebra<EffectiveNumber> EffectiveAlgebra;
typedef ElementaryAlgebra<ApproximateNumber> ApproximateElementaryAlgebra;
typedef ElementaryAlgebra<ValidatedNumber> ValidatedElementaryAlgebra;
typedef ElementaryAlgebra<EffectiveNumber> EffectiveElementaryAlgebra;

template<class F, class P, class SIG> class FunctionMixin { };
template<class F, class P, class... ARGS> class ScalarFunctionMixin;
template<class F, class P, class... ARGS> class VectorFunctionMixin;
template<class F, class P> using ScalarUnivariateFunctionMixin = ScalarFunctionMixin<F,P,RealScalar>;
template<class F, class P> using VectorUnivariateFunctionMixin = VectorFunctionMixin<F,P,RealScalar>;
template<class F, class P> using ScalarMultivariateFunctionMixin = ScalarFunctionMixin<F,P,RealVector>;
template<class F, class P> using VectorMultivariateFunctionMixin = VectorFunctionMixin<F,P,RealVector>;

template<class T> T* heap_copy(const T& t) { return new T(t); }
template<class T> T* heap_move(T&& t) { return new T(std::move(t)); }

template<class P, class D> ScalarFunctionInterface<P,D>* heap_copy(ScalarFunction<P,D> const& f) { return f.raw_pointer()->_clone(); }


template<class D> D make_domain(SizeType d);
template<> inline IntervalDomainType make_domain(SizeType d) { assert(d==1u); return IntervalDomainType(-inf,+inf); }
template<> inline BoxDomainType make_domain(SizeType d) { return BoxDomainType(d,IntervalDomainType(-inf,+inf)); }

template<class F, class SIG>
class FunctionMixin<F,Void,SIG>
    : public virtual FunctionInterface<Void,SIG>
{
    using D=typename SignatureTraits<SIG>::DomainType;
    using C=typename SignatureTraits<SIG>::CodomainType;
  public:
    typedef typename FunctionInterface<Void,SIG>::DomainType DomainType;
    typedef typename FunctionInterface<Void,SIG>::CodomainType CodomainType;
    typedef typename FunctionInterface<Void,SIG>::ArgumentSizeType ArgumentSizeType;
    typedef typename FunctionInterface<Void,SIG>::ResultSizeType ResultSizeType;
  protected:
    FunctionMixin() { }
    template<class X> ElementType<C,X> _base_evaluate(const ElementType<D,X>& x) const;
  public:
    virtual ArgumentSizeType argument_size() const override = 0;
    virtual ResultSizeType result_size() const override = 0;

    virtual OutputStream& _write(OutputStream& os) const override = 0;
    virtual OutputStream& repr(OutputStream& os) const override { return this->_write(os); }
};


template<class F, class SIG>
class FunctionMixin<F,ApproximateTag,SIG>
    : public virtual FunctionInterface<ApproximateTag,SIG>
    , public FunctionMixin<F,Void,SIG>
{
    using D=typename SignatureTraits<SIG>::DomainType;
    using C=typename SignatureTraits<SIG>::CodomainType;
    template<class X> using Argument = typename ElementTraits<D>::template Type<X>;
    template<class X> using Result = typename ElementTraits<C>::template Type<X>;
  public:
    virtual FunctionInterface<ApproximateTag,SIG>* _clone() const override;
    virtual Result<FloatDPApproximation> _evaluate(const Argument<FloatDPApproximation>& x) const override;
    virtual Result<FloatMPApproximation> _evaluate(const Argument<FloatMPApproximation>& x) const override;
    virtual Result<Differential<FloatDPApproximation>> _evaluate(const Argument<Differential<FloatDPApproximation>>& x) const override;
    virtual Result<Differential<FloatMPApproximation>> _evaluate(const Argument<Differential<FloatMPApproximation>>& x) const override;
    virtual Result<TaylorModel<ApproximateTag,FloatDP>> _evaluate(const Argument<TaylorModel<ApproximateTag,FloatDP>>& x) const override;
    virtual Result<TaylorModel<ApproximateTag,FloatMP>> _evaluate(const Argument<TaylorModel<ApproximateTag,FloatMP>>& x) const override;
    virtual Result<Formula<ApproximateNumber>> _evaluate(const Argument<Formula<ApproximateNumber>>& x) const override;
    virtual Result<ElementaryAlgebra<ApproximateNumber>> _evaluate(const Argument<ElementaryAlgebra<ApproximateNumber>>& x) const override;
};

// A wrapper for classes with non-static _compute and _compute_approx methods
template<class F, class SIG>
class FunctionMixin<F,ValidatedTag,SIG>
    : public virtual FunctionInterface<ValidatedTag,SIG>
    , public FunctionMixin<F,ApproximateTag,SIG>
{
    using D=typename SignatureTraits<SIG>::DomainType;
    using C=typename SignatureTraits<SIG>::CodomainType;
    template<class X> using Argument = typename ElementTraits<D>::template Type<X>;
    template<class X> using Result = typename ElementTraits<C>::template Type<X>;
  public:
    using FunctionMixin<F,ApproximateTag,SIG>::_evaluate;
    virtual FunctionInterface<ValidatedTag,SIG>* _clone() const override;
    virtual Result<FloatDPBounds> _evaluate(const Argument<FloatDPBounds>& x) const override;
    virtual Result<FloatMPBounds> _evaluate(const Argument<FloatMPBounds>& x) const override;
    virtual Result<Differential<FloatDPBounds>> _evaluate(const Argument<Differential<FloatDPBounds>>& x) const override;
    virtual Result<Differential<FloatMPBounds>> _evaluate(const Argument<Differential<FloatMPBounds>>& x) const override;
    virtual Result<TaylorModel<ValidatedTag,FloatDP>> _evaluate(const Argument<TaylorModel<ValidatedTag,FloatDP>>& x) const override;
    virtual Result<TaylorModel<ValidatedTag,FloatMP>> _evaluate(const Argument<TaylorModel<ValidatedTag,FloatMP>>& x) const override;
    virtual Result<TaylorModel<ValidatedTag,FloatDPUpperInterval>> _evaluate(const Argument<TaylorModel<ValidatedTag,FloatDPUpperInterval>>& x) const override;
    virtual Result<TaylorModel<ValidatedTag,FloatMPUpperInterval>> _evaluate(const Argument<TaylorModel<ValidatedTag,FloatMPUpperInterval>>& x) const override;
    virtual Result<Formula<ValidatedNumber>> _evaluate(const Argument<Formula<ValidatedNumber>>& x) const override;
    virtual Result<ElementaryAlgebra<ValidatedNumber>> _evaluate(const Argument<ElementaryAlgebra<ValidatedNumber>>& x) const override;

    virtual Result<ValidatedScalarMultivariateFunction> _evaluate(const Argument<ValidatedScalarMultivariateFunction>& x) const override;

};

// A wrapper for classes with non-static _compute and _compute_approx methods
template<class F, class SIG>
class FunctionMixin<F,EffectiveTag,SIG>
    : public virtual FunctionInterface<EffectiveTag,SIG>
    , public FunctionMixin<F,ValidatedTag,SIG>
{
    using D=typename SignatureTraits<SIG>::DomainType;
    using C=typename SignatureTraits<SIG>::CodomainType;
    template<class X> using Argument = typename ElementTraits<D>::template Type<X>;
    template<class X> using Result = typename ElementTraits<C>::template Type<X>;
  public:
    using FunctionMixin<F,ValidatedTag,SIG>::_evaluate;
    virtual FunctionInterface<EffectiveTag,SIG>* _clone() const override;
    virtual Result<Real> _evaluate(const Argument<Real>& x) const override;
    virtual Result<ElementaryAlgebra<Real>> _evaluate(const Argument<ElementaryAlgebra<Real>>& x) const override;
    virtual Result<Formula<Real>> _evaluate(const Argument<Formula<Real>>& x) const override;
    virtual Result<ElementaryAlgebra<EffectiveNumber>> _evaluate(const Argument<ElementaryAlgebra<EffectiveNumber>>& x) const override;
    virtual Result<Formula<EffectiveNumber>> _evaluate(const Argument<Formula<EffectiveNumber>>& x) const override;
};

template<class F, class P, class... ARGS> class ScalarFunctionMixin
    : public FunctionMixin<F,P,RealScalar(ARGS...)> { };

template<class F, class P, class... ARGS> class VectorFunctionMixin
    : public FunctionMixin<F,P,RealVector(ARGS...)>
    , public virtual VectorOfFunctionInterface<P,ARGS...>
{
    virtual ScalarFunctionInterface<P,ARGS...>* _get(SizeType i) const override {
        auto fi=static_cast<F const&>(*this)[i]; return heap_copy(fi); }
};


template<class F,class SIG> FunctionInterface<ApproximateTag,SIG>* FunctionMixin<F,ApproximateTag,SIG>::_clone() const {
    return new F(static_cast<const F&>(*this)); }
template<class F,class SIG> FunctionInterface<ValidatedTag,SIG>* FunctionMixin<F,ValidatedTag,SIG>::_clone() const {
    return new F(static_cast<const F&>(*this)); }
template<class F,class SIG> FunctionInterface<EffectiveTag,SIG>* FunctionMixin<F,EffectiveTag,SIG>::_clone() const {
    return new F(static_cast<const F&>(*this)); }



} // namespace Ariadne

#endif // ARIADNE_FUNCTION_TEMPLATE_HPP
