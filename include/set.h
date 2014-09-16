/***************************************************************************
 *            set.h
 *
 *  Copyright 2013  Pieter Collins
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

/*! \file set.h
 *  \brief Handles for open, closed, overt and compact subsets of Euclidean space.
 */

#ifndef ARIADNE_SET_H
#define ARIADNE_SET_H

#include <iosfwd>

#include "tribool.h"
#include "box.h"
#include "handle.h"
#include "set_interface.h"

namespace Ariadne {

//! \ingroup GeometryModule SetSubModule
//! \brief Base handle class for sets.
class SetBase
    : public Handle<SetInterfaceBase>
{
  public:
    template<class ...Args> SetBase(Args&&... args) : Handle<SetInterfaceBase>(std::forward<Args>(args)...) { }
    inline uint dimension() const { return this->reference().dimension(); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for bounded sets.
class BoundedSet
    : public Handle<BoundedSetInterface>
{
  public:
    template<class ...Args> BoundedSet(Args&&... args) : Handle<BoundedSetInterface>(std::forward<Args>(args)...) { }
    inline uint dimension() const { return this->reference().dimension(); }
    inline Tribool inside(const Box& bx) const { return this->reference().inside(bx); }
    inline Box bounding_box() const { return this->reference().bounding_box(); }
};


//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for overt sets, for which intersection with an open box is verifiable.
class OvertSet
    : public Handle<OvertSetInterface>
{
  public:
    template<typename ...Args> OvertSet(Args&&... args) : Handle<OvertSetInterface>(std::forward<Args>(args)...) { }
    inline uint dimension() const { return this->reference().dimension(); }
    inline Tribool overlaps(const Box& bx) const { return this->reference().overlaps(bx); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for open sets.
class OpenSet
    : public Handle<OpenSetInterface>
{
  public:
    template<typename ...Args> OpenSet(Args&&... args) : Handle<OpenSetInterface>(std::forward<Args>(args)...) { }
    inline uint dimension() const { return this->reference().dimension(); }
    inline Tribool overlaps(const Box& bx) const { return this->reference().overlaps(bx); }
    inline Tribool covers(const Box& bx) const { return this->reference().covers(bx); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for closed sets.
class ClosedSet
    : public Handle<ClosedSetInterface>
{
  public:
    template<typename ...Args> ClosedSet(Args&&... args) : Handle<ClosedSetInterface>(std::forward<Args>(args)...) { }
    inline uint dimension() const { return this->reference().dimension(); }
    inline Tribool separated(const Box& bx) const { return this->reference().separated(bx); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for compact (closed and bounded) sets.
class CompactSet
    : public Handle<CompactSetInterface>
{
  public:
    template<typename ...Args> CompactSet(Args&&... args) : Handle<CompactSetInterface>(std::forward<Args>(args)...) { }
    inline uint dimension() const { return this->reference().dimension(); }
    inline Tribool separated(const Box& bx) const { return this->reference().separated(bx); }
    inline Tribool inside(const Box& bx) const { return this->reference().inside(bx); }
    inline Box bounding_box() const { return this->reference().bounding_box(); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for regular sets, whose closure is the closure of the interior, and whose interior is the interior of the closure.
class RegularSet
    : public Handle<RegularSetInterface>
{
  public:
    template<typename ...Args> RegularSet(Args&&... args) : Handle<RegularSetInterface>(std::forward<Args>(args)...) { }
    inline uint dimension() const { return this->reference().dimension(); }
    inline Tribool overlaps(const Box& bx) const { return this->reference().overlaps(bx); }
    inline Tribool covers(const Box& bx) const { return this->reference().covers(bx); }
    inline Tribool separated(const Box& bx) const { return this->reference().separated(bx); }
};



//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for located (overt and compact) sets.
class LocatedSet
    : public Handle<LocatedSetInterface>
{
  public:
    template<typename ...Args> LocatedSet(Args&&... args) : Handle<LocatedSetInterface>(std::forward<Args>(args)...) { }
    inline uint dimension() const { return this->reference().dimension(); }
    inline Tribool overlaps(const Box& bx) const { return this->reference().overlaps(bx); }
    inline Tribool separated(const Box& bx) const { return this->reference().separated(bx); }
    inline Tribool inside(const Box& bx) const { return this->reference().inside(bx); }
    inline Box bounding_box() const { return this->reference().bounding_box(); }
};

//! \ingroup GeometryModule SetSubModule
//! \brief Handle class for bounded regular sets.
class RegularLocatedSet
    : public Handle<SetInterface>
{
  public:
    template<typename ...Args> RegularLocatedSet(Args&&... args) : Handle<SetInterface>(std::forward<Args>(args)...) { }
    inline uint dimension() const { return this->reference().dimension(); }
    inline Tribool overlaps(const Box& bx) const { return this->reference().overlaps(bx); }
    inline Tribool covers(const Box& bx) const { return this->reference().covers(bx); }
    inline Tribool separated(const Box& bx) const { return this->reference().separated(bx); }
    inline Tribool inside(const Box& bx) const { return this->reference().inside(bx); }
    inline Box bounding_box() const { return this->reference().bounding_box(); }
};



} // namespace Ariadne


#endif // ARIADNE_SET_INTERFACE