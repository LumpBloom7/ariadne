/***************************************************************************
 *            domain.h
 *
 *  Copyright 2008-13  Alberto Casagrande, Pieter Collins
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

/*! \file domain.h
 *  \brief Interval and box domains for functions.
 */

#ifndef ARIADNE_DOMAIN_H
#define ARIADNE_DOMAIN_H

#include "geometry/interval.h"
#include "geometry/box.h"

namespace Ariadne {

using IntervalDomainType = ExactIntervalType;
using BoxDomainType = ExactBoxType;

class RealDomain {
  public:
    constexpr RealDomain() { }
    constexpr SizeOne dimension() const { return SizeOne(); }
    operator IntervalDomainType() const { return IntervalDomainType(-inf,+inf); }
    friend RealDomain intersection(RealDomain const& dom1, RealDomain const& dom2) { return RealDomain(); }
    friend Bool operator==(RealDomain const& dom1, RealDomain const& dom2) { return true; }
    friend OutputStream& operator<<(OutputStream& os, RealDomain const& dom) { return os << "R"; }
};

class EuclideanDomain {
    SizeType _dim;
  public:
    constexpr EuclideanDomain(SizeType dim) : _dim(dim) { }
    constexpr EuclideanDomain(SizeType dim, RealDomain) : _dim(dim) { }
    constexpr SizeType dimension() const { return this->_dim; }
    constexpr RealDomain operator[](SizeType ind) { return RealDomain(); }
    operator BoxDomainType() const { return BoxDomainType(this->dimension(),IntervalDomainType(RealDomain())); }
    friend EuclideanDomain intersection(EuclideanDomain const& dom1, EuclideanDomain const& dom2) { assert(dom1==dom2); return dom1; }
    friend Bool operator==(EuclideanDomain const& dom1, EuclideanDomain const& dom2) { return dom1.dimension() == dom2.dimension(); }
    friend OutputStream& operator<<(OutputStream& os, EuclideanDomain const& dom) { return os << "R" << dom.dimension(); }
};


class UnitInterval {
  public:
    constexpr UnitInterval() { }
    constexpr SizeOne dimension() const { return SizeOne(); }
    operator IntervalDomainType() const { return IntervalDomainType(-1,+1); }
    friend UnitInterval intersection(UnitInterval const& dom1, UnitInterval const& dom2) { return UnitInterval(); }
    friend Bool operator==(UnitInterval const& dom1, UnitInterval const& dom2) { return true; }
    friend OutputStream& operator<<(OutputStream& os, UnitInterval const& dom) { return os << "[-1:+1]"; }
};
class UnitBox {
    SizeType _dim;
  public:
    constexpr UnitBox(SizeType dim) : _dim(dim) { }
    constexpr UnitBox(SizeType dim, UnitInterval) : _dim(dim) { }
    constexpr SizeType dimension() const { return this->_dim; }
    constexpr UnitInterval operator[](SizeType ind) { return UnitInterval(); }
    operator BoxDomainType() const { return BoxDomainType(this->dimension(),IntervalDomainType(UnitInterval())); }
    friend UnitBox intersection(UnitBox const& dom1, UnitBox const& dom2) { assert(dom1==dom2); return dom1; }
    friend Bool operator==(UnitBox const& dom1, UnitBox const& dom2) { return dom1.dimension() == dom2.dimension(); }
    friend OutputStream& operator<<(OutputStream& os, UnitBox const& dom) { return os << "[-1:+1]^" << dom.dimension(); }
};

template<class S> struct ElementTraits;
template<class S, class X> using ElementType = typename ElementTraits<S>::template Type<X>;

template<class UB> struct ElementTraits<Interval<UB>> { template<class X> using Type=Scalar<X>; };
template<class IVL> struct ElementTraits<Box<IVL>> { template<class X> using Type=Vector<X>; };
template<> struct ElementTraits<RealDomain> { template<class X> using Type=Scalar<X>; };
template<> struct ElementTraits<EuclideanDomain> { template<class X> using Type=Vector<X>; };
template<> struct ElementTraits<UnitInterval> { template<class X> using Type=Scalar<X>; };
template<> struct ElementTraits<UnitBox> { template<class X> using Type=Vector<X>; };


} // namespace Ariadne


#endif /* ARIADNE_DOMAIN_H */
