/***************************************************************************
 *            hybrid_set.decl.h
 *
 *  Copyright 2008-17  Pieter Collins
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

/*! \file hybrid_set.decl.h
 *  \brief Forward declarations of sets in hybrid spaces.
 */

#ifndef ARIADNE_HYBRID_SET_DECL_H
#define ARIADNE_HYBRID_SET_DECL_H


namespace Ariadne {

template<class UB> class VariableInterval;
using FloatVariableInterval = VariableInterval<Float64>;


template<class UB> class HybridVariableInterval;
template<class IVL> class HybridVariablesBox;

template<class IVL> class HybridVariablesBox;
template<class IVL> using HybridBoxSet = HybridVariablesBox<IVL>;

using HybridRealBoxSet = HybridVariablesBox<RealInterval>;

class HybridConstraintSet;
class HybridBoundedConstraintSet;
using HybridRealBoundedConstraintSet = HybridBoundedConstraintSet;

using HybridSet = HybridBoundedConstraintSet;

template<class X> class HybridPoint;
using HybridExactPoint = HybridPoint<Float64Value>;


template<class IVL> class HybridBox;
using HybridRealBox = HybridBox<RealInterval>;
using HybridExactBox = HybridBox<ExactIntervalType>;
using HybridUpperBox = HybridBox<ExactIntervalType>;
//using HybridUpperBox = HybridBox<UpperIntervalType>;

using HybridExactBoxType = HybridExactBox;
using HybridUpperBoxType = HybridUpperBox;

template<class IVL> class HybridBoxes;
//using HybridUpperBoxes = HybridBoxes<UpperIntervalType>;
using HybridExactBoxes = HybridBoxes<ExactIntervalType>;
using HybridUpperBoxes = HybridBoxes<ExactIntervalType>;

using HybridExactBoxesType = HybridExactBoxes;
using HybridUpperBoxesType = HybridUpperBoxes;

template<class EBS> class HybridListSet;
class HybridGridTreeSet;

} // namespace Ariadne

#endif // ARIADNE_HYBRID_SET_DECL_H
