/***************************************************************************
 *            geometry/point.hpp
 *
 *  Copyright  2008-20  Alberto Casagrande, Pieter Collins
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

/*! \file geometry/point.hpp
 *  \brief Points in Euclidean space.
 */

#ifndef ARIADNE_POINT_HPP
#define ARIADNE_POINT_HPP

#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "symbolic/identifier.hpp"

#include "output/graphics_interface.hpp"

namespace Ariadne {

template<class X> class Point;
template<class X> class LabelledPoint;

//!@{
//! \relates Point
//! \name Type synonyms
using DyadicPoint = Point<Dyadic>; //!< <p/>
using RationalPoint = Point<Rational>; //!< <p/>
using RealPoint = Point<Real>; //!< <p/>

using LabelledRealPoint = LabelledPoint<Real>;

template<class T> class Space;
using RealSpace = Space<Real>;
template<class T> class Variable;
using RealVariable = Variable<Real>;
template<class V,class E> class Assignment;

template<class F> using ExactPoint = Point<Value<F>>;
template<class F> using ValidatedPoint = Point<Bounds<F>>;
template<class F> using ApproximatePoint = Point<Approximation<F>>;

using FloatDPValuePoint = Point<FloatDPValue>;
using FloatDPBoundsPoint = Point<FloatDPBounds>;
using FloatDPApproximationPoint = Point<FloatDPApproximation>;

typedef ExactPoint<FloatDP> ExactPointType;

template<class IVL> class LabelledBox;
using LabelledExactBoxType = LabelledBox<ExactIntervalType>;
//!@}

//! A point in Euclidean space.
template<class X>
class Point
    : public Vector<X>
    , public DrawableInterface
{
  public:
    typedef X RealType;
    //! A point in zero dimensions
    explicit Point() : Vector<RealType>() { }
    //! The origin in \a n dimensions.
    template<class... PRS> requires Constructible<X,Nat,PRS...>
    explicit Point(SizeType n, PRS... prs) : Vector<RealType>(n,X(0u,prs...)) { }
    Point(const Vector<RealType>& v) : Vector<RealType>(v) { }
    template<class T> requires Convertible<T,X> Point(const Point<T>& pt) : Vector<RealType>(pt.vector()) { }
    template<class Y, class... PRS> requires Constructible<X,Y,PRS...>
    Point(const Point<Y>& pt, PRS... prs) : Vector<RealType>(pt.vector(),prs...) { }
    //! Construct from an initializer list of floating-point values.
    template<class T> requires Convertible<T,X> Point(SizeType n, const T& t) : Vector<RealType>(n,RealType(t)) { }
    //! Construct from an initializer list of floating-point values.
    explicit Point(InitializerList<X> lst);
    //! Construct from a size and an element generator
    template<class G> requires InvocableReturning<X,G,SizeType>
        Point(SizeType n, G const& g) : Vector<X>(n,g) { }
    //! Construct from an initializer list of floating-point values.
    template<class... PRS> requires Constructible<X,ExactDouble,PRS...>
    explicit Point(InitializerList<ExactDouble> lst, PRS... prs);
    //! The origin in \a n dimensions.
    template<class... PRS> requires Constructible<X,Int,PRS...>
    static Point origin(SizeType n, PRS... prs) { return Point(n,RealType(0,prs...)); }
    //! A dynamically-allocated copy.
    virtual Point<X>* clone() const;
    //! The dimension of the point.
    DimensionType dimension() const { return this->size(); }
    //! An explicit cast to a float vector. Useful to prevent ambiguous function overloads.
    const Vector<RealType>& vector() const { return *this; }

    Vector<RealType> centre() const { return *this; }

    friend Point<X> product(Point<X> const& pt1, Point<X> const& pt2) {
        return Point<X>(join(pt1.vector(),pt2.vector())); }

    //! Write to an output stream.
    virtual OutputStream& _write(OutputStream& os) const {
        return os << static_cast<const Vector<RealType>&>(*this); }

    virtual Void draw(CanvasInterface& c, const Projection2d& p) const;
    virtual FloatDPUpperBox bounding_box() const;
};

template<class X>
class LabelledPoint : public LabelledDrawableInterface, public Point<X> {
  public:

    LabelledPoint(Point<X> const& pt, RealSpace const& state_space);
    LabelledPoint(List<Assignment<RealVariable,X>> const& val);
    LabelledPoint(InitializerList<Assignment<RealVariable,X>> const& val);

    LabelledPoint* clone() const override;

    RealSpace state_space() const;

    using Point<X>::draw;
    virtual Void draw(CanvasInterface&, const Variables2d&) const override;

private:
    List<Identifier> _state_variables;
};

template<class X> inline X distance(Point<X> const& pt1, Point<X> const& pt2) {
    ARIADNE_PRECONDITION(pt1.dimension() == pt2.dimension());
    ARIADNE_PRECONDITION(pt1.dimension() > 0);
    X result = sqr(pt1.at(0)-pt2.at(0));
    auto v1=pt1.vector();
    auto v2=pt2.vector();
    for (SizeType i=1; i<pt1.dimension(); ++i) {
        result += sqr(pt1.at(i)-pt2.at(i));
    }
    return sqrt(result);
}

template<class X> Point(Vector<X>) -> Point<X>;

template<class X> template<class... PRS> requires Constructible<X,ExactDouble,PRS...>
Point<X>::Point(InitializerList<ExactDouble> lst, PRS... prs) : Vector<X>(lst,prs...) { }

} // namespace Ariadne

#endif // ARIADNE_POINT_HPP
