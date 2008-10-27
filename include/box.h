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

#include "numeric.h"
#include "vector.h"

#include "set_interface.h"

namespace Ariadne {

class Box
  : public SetInterface,
    public Vector<Interval>
{
 public:
  typedef Float real_type;
  Box() : Vector<Interval>() { }
  template<class T> Box(const T& t) : Vector<Interval>(t) { }
  template<class T1, class T2> Box(const T1& t1, const T2& t2) : Vector<Interval>(t1,t2) { }
  static Box unit_box(uint n) { return Box(n,Interval(-1,1)); }
  static Box upper_quadrant(uint n) { return Box(n,Interval(0,inf())); }
  Vector<Float> centre() const { return midpoint(*this); }
  Float radius() const { 
    Float dmax=0; for(uint i=0; i!=this->size(); ++i) { dmax=max(dmax,(*this)[i].width()); } return up(dmax/2); }

  virtual Box* clone() const { return new Box(*this); }
  virtual uint dimension() const { return this->size(); }
  virtual tribool disjoint(const Vector<Interval>& other) const { 
    return Ariadne::disjoint(*this,other); }
  virtual tribool intersects(const Vector<Interval>& other) const { 
    return !Ariadne::disjoint(*this,other); }
  virtual tribool superset(const Vector<Interval>& other) const { 
    return Ariadne::subset(other,*this); }
  virtual tribool subset(const Vector<Interval>& other) const { 
    return Ariadne::subset(*this,other); }
  virtual Vector<Interval> bounding_box() const { 
    return *this; }
  virtual std::ostream& write(std::ostream& os) const {
    return os << *static_cast<const Vector<Interval>*>(this); }
};

Box make_box(const std::string&);

} // namespace Ariadne

#endif /* ARIADNE_BOX_H */
