/***************************************************************************
 *            polynomial.hpp
 *
 *  Copyright 2008-17  Pieter Collins
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

/*! \file polynomial.hpp
 *  \brief Base class for polynomial rings.
 */

#ifndef ARIADNE_POLYNOMIAL_HPP
#define ARIADNE_POLYNOMIAL_HPP

#include <cassert>
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <algorithm>

#include "../algebra/multi_index.hpp"
#include "../algebra/expansion.hpp"
#include "../algebra/operations.hpp"
#include "../algebra/differential.hpp"


namespace Ariadne {

template<class T> class Array;
template<class X> class Algebra;

//! \brief A monomial with index \a I and coefficients of some type \a X.
template<class X>
class MultivariateMonomial
    : public ExpansionValue<MultiIndex,X>
{
    typedef MultiIndex I;
  public:
    MultivariateMonomial(const MultiIndex& a, const X& x) : ExpansionValue<I,X>(a,x) { }
    MultivariateMonomial(const ExpansionValue<I,X>& v) : ExpansionValue<I,X>(v) { }
};

//! \ingroup FunctionModule
//! \brief A polynomial with coefficients of some type \a X.
template<class X>
class MultivariatePolynomial
    : public DispatchAlgebraOperations<MultivariatePolynomial<X>,X>
{
    template<class XX> friend class MultivariatePolynomial;
    friend struct AlgebraOperations<MultivariatePolynomial<X>,X>;
  public:
    typedef typename Expansion<MultiIndex,X>::ValueType ValueType;
    typedef typename Expansion<MultiIndex,X>::Reference Reference;
    typedef typename Expansion<MultiIndex,X>::ConstReference ConstReference;
    typedef typename Expansion<MultiIndex,X>::Iterator Iterator;
    typedef typename Expansion<MultiIndex,X>::ConstIterator ConstIterator;

    typedef typename Expansion<MultiIndex,X>::IndexReference IndexReference;
    typedef typename Expansion<MultiIndex,X>::IndexConstReference IndexConstReference;
    typedef typename Expansion<MultiIndex,X>::CoefficientReference CoefficientReference;
    typedef typename Expansion<MultiIndex,X>::CoefficientConstReference CoefficientConstReference;

    typedef typename X::Paradigm Paradigm;
    typedef typename X::NumericType NumericType;
    typedef MultivariatePolynomial<X> SelfType;
    typedef ReverseLexicographicIndexLess ComparisonType;
    typedef ReverseLexicographicLess IndexComparisonType;
  public:
    //@{
    //! \name Constructors

    //! \brief The zero polynomial in \a as variables.
    explicit MultivariatePolynomial(SizeType as=0u);
    //! \brief Copy/conversion constructor.
    template<class XX> MultivariatePolynomial(const MultivariatePolynomial<XX>& p);
    //! \brief Copy/conversion constructor.
    template<class XX> explicit MultivariatePolynomial(const Expansion<MultiIndex,XX>& e);
    //! \brief A sparse polynomial with coefficients given by an initializer list of indices and coefficients.
    MultivariatePolynomial(InitializerList<Pair<InitializerList<DegreeType>,X>> lst);
    //@}

    //! \brief Create the null polynomial in the same number of variables.
    MultivariatePolynomial<X> create_zero() const;

    //! \brief Create a constant polynomial in \a as variables with value \a c.
    static MultivariatePolynomial<X> constant(SizeType as, const X& c);
    //! \brief Create a polynomial in \a as variables which returns the value of the \a j<sup>th</sup> variable.
    static MultivariatePolynomial<X> coordinate(SizeType as, SizeType j);
    static MultivariatePolynomial<X> variable(SizeType as, SizeType j);
    //! \brief Create an Array of polynomials in \a as variables,
    //! the i<sup>th</sup> of  which returns the value of the i<sup>th</sup> variable.
    static Vector<MultivariatePolynomial<X>> coordinates(SizeType as);
    static Vector<MultivariatePolynomial<X>> variables(SizeType as);

    //! \brief Set equal to a constant.
    MultivariatePolynomial<X>& operator=(const X& x);
    //@{
    //! \name Comparisons

    //! \brief Equality operator.
    template<class XX> EqualityType<X,XX> operator==(const MultivariatePolynomial<XX>& p) const;
    //! \brief Inequality operator.
    template<class XX> InequalityType<X,XX> operator!=(const MultivariatePolynomial<XX>& p) const;
    //@}

    //@{
    //! \name Data access

    //! \brief The number of variables in the argument of the polynomial.
    SizeType argument_size() const;
    //! \brief The number of structural nonzero terms.
    SizeType number_of_terms() const;
    //! \brief The order of the highest term.
    DegreeType degree() const;
    //! \brief The value of the polynomial at zero.
    const X& value() const;
    //! \brief A reference to the coefficient of the term in \f$x^{a_1}\cdots x^{a_n}\f$.
    X& operator[](const MultiIndex& a);
    //! \brief A constant referent to the coefficient of the term in \f$x^{a_1}\cdots x^{a_n}\f$.
    const X& operator[](const MultiIndex& a) const;
    //! \brief A constant reference to the raw data expansion.
    const Expansion<MultiIndex,X>& expansion() const;
    //! \brief A reference to the raw data expansion.
    Expansion<MultiIndex,X>& expansion();
    //@}

    //@{
    //! \name Iterators

    //! \brief An Iterator to the beginning of the list of terms.
    Iterator begin();
    //! \brief An Iterator to the end of the list of terms..
    Iterator end();
    //! \brief An Iterator to the term in \f$x^a\f$. Returns \c end() if there is no term in \a a.
    Iterator find(const MultiIndex& a);
    //! \brief A constant Iterator to the beginning of the list of terms.
    ConstIterator begin() const;
    //! \brief A constant Iterator to the end of the list of terms.
    ConstIterator end() const;
    //! \brief An Iterator to the term in \f$x^a\f$. Returns \c end() if there is no term in \a a.
    ConstIterator find(const MultiIndex& a) const;
    //@}


    //@{
    //! \name Modifying operations

    //! \brief Insert the term \f$c x^{a_1}\f$ into a sorted list of terms.
    Void insert(const MultiIndex& a, const X& c);
    //! \brief Reserve space for a total of \a n terms.
    Void reserve(SizeType n);
    //! \brief Remove the term pointed to by \a iter. May be expensive if the term is near the beginning of the list of terms.
    Void erase(Iterator iter);
    //! \brief Set the polynomial to zero.
    Void clear();
    //! \brief Remove all zero terms from the expansion, and order the expansion reverse lexicographically by term.
    Void cleanup();
    //@}

    //@{
    //! \name Evaluation

    //! Evaluate on a vector of algebra elements.
    template<class A> A operator() (Vector<A> const&) const;
    //@}

    //@{
    //! \name Modifying operators

    //! \brief Truncate to degree \a d.
    MultivariatePolynomial<X>& truncate(DegreeType d);
    //! \brief Differentiate with respect to the \a j<sup>th</sup> variable.
    MultivariatePolynomial<X>& differentiate(SizeType j);
    //! \brief Antidifferentiate (integrate) with respect to the \a j<sup>th</sup> variable.
    MultivariatePolynomial<X>& antidifferentiate(SizeType j);
    //@}

    //@{
    //! \name Related operations
    friend MultivariatePolynomial<X>& operator*=(MultivariatePolynomial<X>& p, const MultivariateMonomial<X>& m) { return MultivariatePolynomial<X>::_imul(p,m); }

    template<class XX, class A> friend A evaluate(const MultivariatePolynomial<XX>& p, const Vector<A>& v);
    template<class XX> friend MultivariatePolynomial<XX> compose(const MultivariatePolynomial<XX>& p, const Vector<MultivariatePolynomial<XX>>& q);
    template<class XX> friend MultivariatePolynomial<XX> derivative(MultivariatePolynomial<XX> dx, SizeType k);
    template<class XX> friend MultivariatePolynomial<XX> antiderivative(MultivariatePolynomial<XX> dx, SizeType k);
    template<class XX> friend MultivariatePolynomial<XX> truncate(MultivariatePolynomial<XX> dx, DegreeType deg);
    //@}

    Void check() const;
    static MultivariatePolynomial<X> _compose(const MultivariatePolynomial<X>& p, const Vector<MultivariatePolynomial<X>>& q);
    static X _evaluate(const MultivariatePolynomial<X>& p, const Vector<X>& vx);
    static Algebra<X> _evaluate(const MultivariatePolynomial<X>& p, const Vector<Algebra<X>>& va);
    static MultivariatePolynomial<X> _partial_evaluate(const MultivariatePolynomial<X>& p, SizeType k, const X& c);
    OutputStream& _write(OutputStream& os) const;
    OutputStream& _write(OutputStream& os, List<String> const& names) const;
  private:
    Void _append(const MultiIndex& a, const X& c);
    Iterator _unique_key();
  private:
    SortedExpansion<MultiIndex,X,ReverseLexicographicIndexLess> _expansion;
  private: // FIXME: Put these concrete-generic operations in proper place
    template<class Y, EnableIf<IsAssignable<X,Y>> =dummy>
        friend MultivariatePolynomial<X> operator+(MultivariatePolynomial<X> p, const Y& c) {
            X xc=p.value(); xc=c; return p+xc; }
    template<class Y, EnableIf<IsAssignable<X,Y>> =dummy>
        friend MultivariatePolynomial<X> operator-(MultivariatePolynomial<X> p, const Y& c) {
            X xc=p.value(); xc=c; return p-xc; }
    template<class Y, EnableIf<IsAssignable<X,Y>> =dummy>
        friend MultivariatePolynomial<X> operator*(const Y& c, MultivariatePolynomial<X> p) {
            X xc=p.value(); xc=c; return xc*p; }
    template<class Y, EnableIf<IsAssignable<X,Y>> =dummy>
        friend MultivariatePolynomial<X> operator*(MultivariatePolynomial<X> p, const Y& c) {
            X xc=p.value(); xc=c; return p*xc; }
    template<class Y, EnableIf<IsAssignable<X,Y>> =dummy>
        friend MultivariatePolynomial<X> operator/(MultivariatePolynomial<X> p, const Y& c) {
            X xc=p.value(); xc=c; return p/xc; }

};

template<class X> struct AlgebraOperations<MultivariatePolynomial<X>> {
  public:
    static MultivariatePolynomial<X> apply(Pos, const MultivariatePolynomial<X>& p);
    static MultivariatePolynomial<X> apply(Neg, const MultivariatePolynomial<X>& p);
    static MultivariatePolynomial<X> apply(Add, const MultivariatePolynomial<X>& p1, const MultivariatePolynomial<X>& p2);
    static MultivariatePolynomial<X> apply(Sub, const MultivariatePolynomial<X>& p1, const MultivariatePolynomial<X>& p2);
    static MultivariatePolynomial<X> apply(Mul, const MultivariatePolynomial<X>& p1, const MultivariatePolynomial<X>& p2);
    static MultivariatePolynomial<X> apply(Add, MultivariatePolynomial<X> p, const X& c);
    static MultivariatePolynomial<X> apply(Mul, MultivariatePolynomial<X> p, const X& c);
    static MultivariatePolynomial<X> apply(Mul, MultivariatePolynomial<X> p, const MultivariateMonomial<X>& m);
    static MultivariatePolynomial<X>& iapply(Add, MultivariatePolynomial<X>& p, const X& c);
    static MultivariatePolynomial<X>& iapply(Mul, MultivariatePolynomial<X>& p, const X& c);
    static MultivariatePolynomial<X>& iapply(Mul, MultivariatePolynomial<X>& p, const MultivariateMonomial<X>& m);

};


template<class X> template<class XX> MultivariatePolynomial<X>::MultivariatePolynomial(const MultivariatePolynomial<XX>& p)
    : _expansion(p._expansion) { }

template<class X> template<class XX> MultivariatePolynomial<X>::MultivariatePolynomial(const Expansion<MultiIndex,XX>& e)
    : _expansion(e) { this->cleanup(); }

template<class X> template<class XX> EqualityType<X,XX> MultivariatePolynomial<X>::operator==(const MultivariatePolynomial<XX>& p) const {
    const_cast<MultivariatePolynomial<X>*>(this)->cleanup();
    const_cast<MultivariatePolynomial<XX>&>(p).cleanup();
    return this->_expansion==p._expansion;
}

template<class X> template<class XX> InequalityType<X,XX> MultivariatePolynomial<X>::operator!=(const MultivariatePolynomial<XX>& p) const {
    return !(*this==p);
}

template<class X> inline MultivariatePolynomial<X> partial_evaluate(const MultivariatePolynomial<X>& p, SizeType k, const X& c) {
    return MultivariatePolynomial<X>::_partial_evaluate(p,k,c); }

template<class X> inline X evaluate(const MultivariatePolynomial<X>& p, const Vector<X>& v) {
    return MultivariatePolynomial<X>::_evaluate(p,v); }

template<class X> inline MultivariatePolynomial<X> compose(const MultivariatePolynomial<X>& p, const Vector<MultivariatePolynomial<X>>& q) {
    return MultivariatePolynomial<X>::_compose(p,q); }

template<class X> inline Vector<MultivariatePolynomial<X>> compose(const Vector<MultivariatePolynomial<X>>& p, const Vector<MultivariatePolynomial<X>>& q) {
    return MultivariatePolynomial<X>::_compose(p,q); }

template<class X, class A> inline A evaluate(const MultivariatePolynomial<X>& p, const Vector<A>& v) {
    return horner_evaluate(p.expansion(),v); }

template<class X> inline MultivariatePolynomial<X> derivative(MultivariatePolynomial<X> p, SizeType k) {
    p.differentiate(k); return std::move(p); }

template<class X> inline MultivariatePolynomial<X> antiderivative(MultivariatePolynomial<X> p, SizeType k) {
    p.antidifferentiate(k); return std::move(p); }

template<class X> inline MultivariatePolynomial<X> truncate(MultivariatePolynomial<X> p, DegreeType deg) {
    p.truncate(deg); return std::move(p); }

template<class X> inline OutputStream& operator<<(OutputStream& os, const MultivariatePolynomial<X>& p) {
    return p._write(os); }

template<class X> inline Bool compatible(const MultivariatePolynomial<X>& x1, const MultivariatePolynomial<X>& x2) {
    return x1.argument_size()==x2.argument_size(); }

template<class F> struct NamedArgumentRepresentation {
    const F& function; const List<String>& argument_names;
};

template<class F> inline NamedArgumentRepresentation<F> named_argument_repr(const F& function, const List<String>& argument_names) {
    NamedArgumentRepresentation<F> r={function,argument_names}; return r; }

template<class X> inline OutputStream& operator<<(OutputStream& os, const NamedArgumentRepresentation<MultivariatePolynomial<X>>& repr) {
    return repr.function._write(os,repr.argument_names); }

template<class X> inline MultivariatePolynomial<MidpointType<X>> midpoint(const MultivariatePolynomial<X>& p) {
    return MultivariatePolynomial<MidpointType<X>>(midpoint(p.expansion())); }


// Vectorised operations
template<class X, class A> Vector<A> evaluate(const Vector<MultivariatePolynomial<X>>& p, const Vector<A>& v) {
    Vector<A> r(p.size(),v.zero_element());
    for(Nat i=0; i!=p.size(); ++i) { r[i]=evaluate(p[i],v); }
    return r;
}

template<class X> Vector<MultivariatePolynomial<X>> derivative(const Vector<MultivariatePolynomial<X>>& p, SizeType j) {
    Vector<MultivariatePolynomial<X>> r(p.size());
    for(Nat i=0; i!=p.size(); ++i) { r[i]=derivative(p[i],j); }
    return r;
}

template<class X> Vector<MultivariatePolynomial<X>> antiderivative(const Vector<MultivariatePolynomial<X>>& p, SizeType j) {
    Vector<MultivariatePolynomial<X>> r(p.size());
    for(Nat i=0; i!=p.size(); ++i) { r[i]=antiderivative(p[i],j); }
    return r;
}

template<class X> Vector<MultivariatePolynomial<X>> truncate(const Vector<MultivariatePolynomial<X>>& p, DegreeType d) {
    Vector<MultivariatePolynomial<X>> r(p.size());
    for(Nat i=0; i!=p.size(); ++i) { r[i]=truncate(p[i],d); }
    return r;
}

template<class X> Vector<MultivariatePolynomial<MidpointType<X>>> midpoint(const Vector<MultivariatePolynomial<X>>& p) {
    Vector<MultivariatePolynomial<MidpointType<X>>> r(p.size());
    for(Nat i=0; i!=p.size(); ++i) { r[i]=midpoint(p[i]); }
    return r;
}


} // namespace Ariadne

#endif /* ARIADNE_POLYNOMIAL_HPP */
