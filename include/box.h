/***************************************************************************
 *            box.h
 *
 *  Copyright 2008  Alberto Casagrande, Pieter Collins
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

/*! \file box.h
 *  \brief Coordinate-aligned boxes in Euclidean space.
 */

#ifndef ARIADNE_BOX_H
#define ARIADNE_BOX_H

#include "container.h"

#include "numeric.h"
#include "vector.h"

#include "set_interface.h"
#include "graphics_interface.h"
#include "point.h"

namespace Ariadne {

typedef std::ostream OutputStream;
typedef uint Nat;

class IntervalSet;
class BoxSet;

template<class X> class Point;
typedef Point<ExactFloatType> ExactPoint;
class Box;

bool contains(const Vector<Interval>& v1, const Vector<ExactFloatType>& v2);
bool contains(const Vector<Interval>& v1, const Vector<ValidatedFloatType>& v2);
bool element(const Vector<ValidatedFloatType>& v1, const Vector<Interval>& v2);

//! \ingroup GeometryModule ExactSetSubModule
//! \brief An exact interval in \f$\mathbb{R}\f$.
class IntervalSet {
    Real _lower, _upper;
  public:
    IntervalSet() : _lower(-1), _upper(+1) { }
    IntervalSet(const Real& l, const Real& u) : _lower(l), _upper(u) { }
    const Real& lower() const { return _lower; }
    const Real& upper() const { return _upper; }
    const Real midpoint() const { return (_lower+_upper)/2; }
    const Real radius() const { return (_upper-_lower)/2; }
};
inline OutputStream& operator<<(OutputStream& os, const IntervalSet& ivl) {
    return os << "{" << ivl.lower() << ":" << ivl.upper() << "}";
}
inline Interval under_approximation(const IntervalSet& rivl) {
    return Interval(Interval(rivl.lower()).upper(),Interval(rivl.upper()).lower());
}
inline Interval over_approximation(const IntervalSet& rivl) {
    return Interval(Interval(rivl.lower()).lower(),Interval(rivl.upper()).upper());
}
inline Interval approximation(const IntervalSet& rivl) {
    return Interval(RawFloatType(ApproximateNumberType(rivl.lower())),RawFloatType(ApproximateNumberType(rivl.upper())));
}


//! \ingroup GeometryModule ExactSetSubModule
//! \brief An exact coordinate-aligned box in \f$\mathbb{R}^n\f$.
class BoxSet {
    Array<IntervalSet> _ary;
  public:
    BoxSet() : _ary() { }
    explicit BoxSet(const Vector<Interval>& iv);
    BoxSet(const List<IntervalSet>& t) : _ary(t.begin(),t.end()) { }
    BoxSet(Nat n, const IntervalSet& ivl) : _ary(n,ivl) { }
    Nat size() const { return _ary.size(); }
    Nat dimension() const { return _ary.size(); }
    IntervalSet const& operator[](Nat i) const { return _ary[i]; }
    IntervalSet& operator[](Nat i) { return _ary[i]; }
    friend OutputStream& operator<<(OutputStream& os, const BoxSet& bx) { return os << bx._ary; }
};
Box under_approximation(const BoxSet& rbx);
Box over_approximation(const BoxSet& rbx);
Box approximation(const BoxSet& rbx);

Box widen(const Box& bx);

//! \ingroup BasicSetSubModule GeometryModule
//! \brief A box in Euclidean space.
class Box
    : public SetInterface,
      public DrawableInterface,
      public Vector<Interval>
{
  public:
    //! Construct a singleton point in zero dimensions.
    Box() : Vector<Interval>() { }
    //! Construct an empty box in \a d dimensions.
    explicit Box(uint d) : Vector<Interval>(d) { }
    //! Construct from an initializer list of pairs of floating-point values
    //! giving lower and upper bounds.
    Box(std::initializer_list<Interval> lst);

    explicit Box(uint d, Interval ivl) : Vector<Interval>(d,ivl) { }
    explicit Box(const Vector<ValidatedFloat>& vec) : Vector<Interval>(vec) { }
    Box(const Vector<Interval>& ivec) : Vector<Interval>(ivec) { }
    Box(const List<Interval>& ilst) : Vector<Interval>(ilst) { }

    //! Construct from a string literal of the form "[a1,b1]x[a2,b2]x...x[ad,bd]".
    explicit Box(const std::string& str);

    //! The unit box \f$[-1,1]^n\f$ in \a n dimensions.
    static Box unit_box(uint n) {
        return Box(n,Interval(-1,1));
    }

    //! The upper quadrant box \f$[0,\infty]^n\f$ in \a n dimensions.
    static Box upper_quadrant(uint n) {
        return Box(n,Interval(0,inf));
    }

    //! An explicit case to an interval vector. Useful to prevent ambiguous function overloads.
    const Vector<Interval>& vector() const { return *this; }

    //! The set of vertices of the box.
    std::vector<ExactPoint> vertices() const;

    //! An approximation to the centre of the box.
    ExactPoint centre() const {
        return ExactPoint(midpoint(static_cast<const Vector<Interval>&>(*this)));
    }

    //! An over-approximation to radius of the box in the supremum norm.
    ErrorType radius() const {
        ErrorType dmax=0;
        for(uint i=0; i!=this->size(); ++i) {
            dmax = max( dmax, (*this)[i].width() );
        }
        return half(dmax);
    }

    //! An approximation to the Lesbegue measure (area, volume) of the box.
    ApproximateNumberType measure() const {
        ApproximateNumberType meas=1;
        for(uint i=0; i!=this->size(); ++i) {
            meas *= (*this)[i].width();
        }
        return meas;
    }

    //! \brief Test if the box is empty.
    bool empty() const {
        return Ariadne::empty(*this);
    }

    //! \brief Test if the box is bounded.
    bool bounded() const {
        for(uint i=0; i!=this->Vector<Interval>::size(); ++i) {
            if(!Ariadne::bounded((*this)[i])) {
                return false;
            }
        }
        return true;
    }

    //! \brief Make a dynamically-allocated copy.
    virtual Box* clone() const {
        return new Box(*this);
    }

    //! \brief Test if the box is a superset of another box.
    bool contains(const ExactPoint& pt) const {
        return Ariadne::contains(this->vector(),pt.vector());
    }

    //! \brief Test if the box is a subset of another box.
    bool subset(const Box& bx) const {
        return Ariadne::subset(this->vector(),bx.vector());
    }

    //! \brief Test if the box is a superset of another box.
    bool superset(const Box& bx) const {
        return Ariadne::subset(bx.vector(),this->vector());
    }

    //! \brief Test if the box intersects another box. Returns \a true even
    //! if the intersection occurs only on the boundary of the boxes.
    //! Use Box::overlaps(const Box& bx) to test robust intersection
    //! of the interiors.
    bool intersects(const Box& bx) const {
        return Ariadne::intersect(this->vector(),bx.vector());
    }

    //! \brief Test if the box intersects another box. Returns \a false even
    //! if the intersection only occurs only on the boundary of the boxes.
    bool disjoint(const Box& bx) const {
        return Ariadne::disjoint(this->vector(),bx.vector());
    }

    //! \brief The dimension of the space the box lies in.
    virtual uint dimension() const {
        return this->size();
    }

    //! \brief Tests if the box is disjoint from another box.
    //! Only returns true if the closures are disjoint.
    virtual tribool separated(const Box& other) const {
        return Ariadne::disjoint(this->vector(), other.vector());
    }

    //! \brief Tests if the box overlaps another box.
    //! Only returns true if the interior of one box intersects the other.
    virtual tribool overlaps(const Box& other) const {
        return Ariadne::overlap(this->vector(), other.vector());
    }

    //! \brief Tests if the box covers another box.
    //! Only returns true if the interior of the box is a superset
    //! of the closure of the other.
    virtual tribool covers(const Box& other) const {
        return Ariadne::inside(other.vector(), this->vector());
    }

    //! \brief Tests if the box covers another box.
    //! Only returns true if the closure of the box is a subset
    //! of the interior of the other.
    virtual tribool inside(const Box& other) const {
        return Ariadne::inside(this->vector(), other.vector());
    }

    //! \brief Returns an enclosing bounding box for the set.
    //! The result is guaranteed to have nonempty interior, and floating-point
    //! boundary coefficients so the centre and radius are exactly computable.
    virtual Box bounding_box() const {
        return Ariadne::widen(*this);
    }

    //! \brief Widens the box by the minimal floating-point increment.
    //! The result is guaranteed to contain the box in its interior.
    Void widen() {
        Box& result=*this;
        for(uint i=0; i!=result.dimension(); ++i) {
            result[i]=Ariadne::widen((*this)[i]);
        }
    }

    //! \brief Narrows the box by the minimal floating-point increment.
    //! The result is guaranteed to contain the box in its interior.
    Void narrow() {
        Box& result=*this;
        for(uint i=0; i!=result.dimension(); ++i) {
            result[i]=Ariadne::narrow((*this)[i]);
        }
    }

    //! \brief Split into two along the largest side.
    std::pair<Box,Box> split() const { return Ariadne::split(*this); }
    //! \brief Split into two along side with index \a i.
    std::pair<Box,Box> split(uint i) const { return Ariadne::split(*this,i); };

    //! \brief Draw on a canvas.
    virtual void draw(CanvasInterface& c, const Projection2d& p) const;

    //! \brief Write to an output stream.
    virtual std::ostream& write(std::ostream& os) const {
        return os << *static_cast<const Vector<Interval>*>(this);
    }

    operator UpperBox & () { return reinterpret_cast<UpperBox&>(static_cast<Vector<Interval>&>(*this)); }
    operator UpperBox const& () const { return reinterpret_cast<UpperBox const&>(static_cast<Vector<Interval>const&>(*this)); }
    operator Vector<UpperInterval> const& () const { return reinterpret_cast<Vector<UpperInterval>const&>(static_cast<Vector<Interval>const&>(*this)); }
};

//! \relates Box \brief The cartesian product of two boxes.
Box product(const Box& bx1, const Box& bx2);
//! \relates Box \brief The smallest box containing the two boxes.
Box hull(const Box& bx1, const Box& bx2);
Box hull(const Box& bx1, const ExactPoint& pt2);
Box hull(const ExactPoint& pt1, const ExactPoint& pt2);
//! \relates Box \brief The intersection of the two boxes.
Box intersection(const Box& bx1, const Box& bx2);
//! \relates Box \brief A box which is wider than the input, and has single-precision values.
Box widen(const Box& bx);
//! \relates Box \brief A box which is narrower than the input, and has single-precision values.
Box narrow(const Box& bx);

Box make_box(const std::string& str);


//! \ingroup BasicSetSubModule GeometryModule
//! \brief A box in Euclidean space.
class UpperBox
    : public Vector<UpperInterval>
{
  public:
    //! Construct a singleton point in zero dimensions.
    UpperBox() : Vector<UpperInterval>() { }
    //! Construct an empty box in \a d dimensions.
    explicit UpperBox(uint d) : Vector<UpperInterval>(d,UpperInterval()) { }
    explicit UpperBox(uint d, UpperInterval ivl) : Vector<UpperInterval>(d,ivl) { }
    //! Construct from an initializer list of pairs of floating-point values
    //! giving lower and upper bounds.
    UpperBox(std::initializer_list<UpperInterval> lst);

    UpperBox(Vector<Interval>const& vec) : Vector<UpperInterval>(vec) { }
    UpperBox(Vector<UpperInterval>const& vec) : Vector<UpperInterval>(vec) { }
    template<class E> UpperBox(const VectorExpression<E>& t) : Vector<UpperInterval>(t) { }

    explicit UpperBox(const Vector<ValidatedFloat>& vec) : Vector<UpperInterval>(vec) { }

    //! The unit box \f$[-1,1]^n\f$ in \a n dimensions.
    static UpperBox unit_box(uint n) {
        return UpperBox(n,UpperInterval(-1.0,1.0));
    }

    //! The upper quadrant box \f$[0,\infty]^n\f$ in \a n dimensions.
    static UpperBox upper_quadrant(uint n) {
        return UpperBox(n,UpperInterval(0,inf));
    }

    //! An explicit case to an interval vector. Useful to prevent ambiguous function overloads.
    const Vector<UpperInterval>& vector() const { return *this; }

    //! An approximation to the centre of the box.
    ApproximatePoint centre() const {
        ApproximatePoint r(this->dimension());
        for(uint i=0; i!=this->dimension(); ++i) { r[i]=midpoint((*this)[i]); }
        return r;
    }

    //! An over-approximation to radius of the box in the supremum norm.
    ErrorType radius() const {
        ErrorType dmax=0;
        for(uint i=0; i!=this->size(); ++i) { dmax = max( dmax, (*this)[i].width() ); }
        return half(dmax);
    }

    //! An approximation to the Lesbegue measure (area, volume) of the box.
    ApproximateNumberType measure() const {
        ApproximateNumberType meas=1;
        for(uint i=0; i!=this->size(); ++i) { meas *= (*this)[i].width(); }
        return meas;
    }

    //! \brief Test if the box is bounded.
    bool bounded() const {
        for(uint i=0; i!=dimension(); ++i) {
            if(!Ariadne::bounded((*this)[i])) { return false; } }
        return true;
    }

    //! \brief Test if the box is a subset of another box.
    bool refines(const UpperBox& bx) const {
        for(uint i=0; i!=dimension(); ++i) {
            if(!Ariadne::refines((*this)[i],bx[i])) { return false; } }
        return true;
    }

    friend bool refines(const UpperBox& bx1, const UpperBox& bx2) {
        return bx1.refines(bx2);
    }

    friend UpperBox intersection(const UpperBox& bx1, const UpperBox& bx2) {
        assert(bx1.dimension()==bx2.dimension());
        UpperBox rbx(bx1.dimension());
        for(uint i=0; i!=rbx.dimension(); ++i) { rbx[i]=intersection(bx1[i],bx2[i]); }
        return rbx;
    }

    //! \brief Test if the box intersects another box. Returns \a false even
    //! if the intersection only occurs only on the boundary of the boxes.
    Tribool disjoint(const UpperBox& bx) const {
        for(uint i=0; i!=this->dimension(); ++i) {
            if(Ariadne::disjoint((*this)[i],bx[i])) { return true; } }
        return indeterminate;
    }

    friend Tribool disjoint(const UpperBox& bx1, const UpperBox& bx2) {
        return bx1.disjoint(bx2);
    }

    //! \brief The dimension of the space the box lies in.
    uint dimension() const {
        return this->size();
    }

    //! \brief Tests if the box is disjoint from another box.
    //! Only returns true if the closures are disjoint.
    tribool empty() const{
        for(uint i=0; i!=this->dimension(); ++i) {
            if((*this)[i].empty()) { return true; } }
        return indeterminate;
    }

    //! \brief Tests if the box is disjoint from another box.
    //! Only returns true if the closures are disjoint.
    tribool separated(const UpperBox& other) const {
        return this->disjoint(other);
    }

    //! \brief Returns an enclosing bounding box for the set.
    //! The result is guaranteed to have nonempty interior, and floating-point
    //! boundary coefficients so the centre and radius are exactly computable.
    UpperBox bounding_box() const {
        UpperBox result=*this; result.widen(); return result;
    }

    //! \brief Widens the box by the minimal floating-point increment.
    //! The result is guaranteed to contain the box in its interior.
    Void widen() {
        UpperBox& result=*this;
        for(uint i=0; i!=result.dimension(); ++i) {
            result[i]=Ariadne::widen((*this)[i]);
        }
    }

    //! \brief Split into two along the largest side.
    std::pair<UpperBox,UpperBox> split() const;
    //! \brief Split into two along side with index \a i.
    std::pair<UpperBox,UpperBox> split(uint i) const;

    //! \brief Write to an output stream.
    std::ostream& write(std::ostream& os) const {
        return os << *static_cast<const Vector<UpperInterval>*>(this);
    }
};

inline Box make_exact_box(UpperBox const& ubx) { return Box(reinterpret_cast<Vector<Interval>const&>(static_cast<Vector<UpperInterval>const&>(ubx))); }
Tribool disjoint(const UpperBox& bx1, const UpperBox& bx2);

} // namespace Ariadne

#endif /* ARIADNE_BOX_H */
