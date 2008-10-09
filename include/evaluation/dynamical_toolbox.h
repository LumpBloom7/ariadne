/***************************************************************************
 *            dynamical_toolbox.h
 *
 *  Copyright  2008  Pieter Collins
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
 
/*! \file dynamical_toolbox.h
 *  \brief Toolbox for dynamical systems. 
 */


#ifndef ARIADNE_DYNAMICAL_TOOLBOX_H
#define ARIADNE_DYNAMICAL_TOOLBOX_H

#include "base/tuple.h"

/* \brief Top-level namespace. */
namespace Ariadne {
 
using std::pair;

template<class T, ushort N> class array;

template<class R> class Interval;
template<class R> class FunctionInterface;
template<class X> class Vector;
template<class R> class Box;
template<class R> class TaylorSet;
class DiscreteEvent;

enum CrossingKind { MISSING, TOUCHING, TRANSVERSE, CROSSING, GRAZING, UNKNOWN };


/*! \brief Tools for analysing dynamical systems based on function models. */
template<class Mdl> 
class DynamicalToolbox
{
  typedef typename Mdl::real_type R;
  typedef typename traits<R>::approximate_arithmetic_type A;
  typedef Interval<R> I;

 public:
  //!
  typedef R RealType;
  //!
  typedef Mdl ModelType;
  typedef TaylorSet<R> SetType;
  //!
  typedef Mdl SetModelType;
  typedef Mdl TimeModelType;
  typedef Mdl MapModelType;
  typedef Mdl FlowModelType;
  typedef Mdl GuardModelType;
  typedef R TimeType;
  //typedef Box<RealType> BoxType;
  typedef Vector<I> BoxType;
  typedef Interval<RealType> IntervalType;
  typedef FunctionInterface<RealType> FunctionType;
 public:
  //! \brief Test if a box satisfies the constraint given by the guard. Returns \a true is all points
  //! in the box satisfy the constraint, \a false if all points do not satisfy the constraint, and 
  //! indeterminate otherwise.
  tribool active(const FunctionType& guard,  const BoxType& box) const;

  //! \brief Test if a set satisfies the constraint given by the guard. Returns \a true is all points
  //! in the set satisfy the constraint, \a false if all points do not satisfy the constraint, and 
  //! indeterminate otherwise.
  tribool active(const FunctionType& guard,  const SetModelType& set_model) const;

  //! \brief Test if a set satisfies the constraint given by the guard model. Returns \a true is all 
  //! points in the set satisfy the constraint, \a false if all points do not satisfy the constraint, 
  //! and indeterminate otherwise.
  tribool active(const GuardModelType& guard_model, 
                 const SetModelType& set_model) const;

  //! \brief Computes an over-approximation to the time interval for which the \a initial_set_model 
  //! touch the set specified by the \a guard model under the \a flow_model. The \a minimum and \a maximum_time 
  //! gives the minimum and maximum time for which the evolution is valid.
  pair<TimeModelType,TimeModelType> 
  touching_time_interval(const FlowModelType& flow_model, 
                         const GuardModelType& guard_model, 
                         const SetModelType& initial_set_model) const;
  
  //! \brief Computes the time at which points in the \a initial_set_model cross the zero-set of the
  //! the \a guard_model under evolution of the \a flow_model, for times between the \a minimum_time and \a maximum_time.
  //! The crossing must be (differentiably) transverse.
  TimeModelType crossing_time(const FlowModelType& flow_model, 
                              const GuardModelType& guard_model, 
                              const SetModelType& initial_set_model) const;

  
  
  tuple<DiscreteEvent,CrossingKind,Mdl,Mdl>
  crossing_data(const FlowModelType& flow_model, 
                const GuardModelType& guard_model, 
                const SetModelType& initial_set_model) const;

  //! \brief Computes the image of the set defined by \a set_model under the \a map.
  SetModelType reset_step(const FunctionType& map, 
                          const SetModelType& set_model) const;
  
  //! \brief Computes the image of the set defined by \a set_model under the approximation of the map 
  //! given by \a map_model.
  SetModelType reset_step(const MapModelType& map_model, 
                          const SetModelType& set_model) const;
  
  //! \brief Computes the points reached by evolution of the \a initial_set_model under the flow
  //! given by \a flow_model. The \a integration_time gives the time all points should be flowed.
  SetModelType integration_step(const FlowModelType& flow_model, 
                                const SetModelType& initial_set_model, 
                                const TimeType& integration_time) const;
  
  //! \brief Computes the points reached by evolution of the \a initial_set_model under the flow
  //! given by \a flow_model. The \a integration_time_model \f$\tau(e)\f$ gives the time the point 
  //! starting at \f$x(e)\f$ should be flowed.
  SetModelType integration_step(const FlowModelType& flow_model, 
                                const SetModelType& initial_set_model, 
                                const TimeModelType& integration_time_model) const;
  
  //! \brief Computes the points reached by evolution of the \a initial_set_model under the flow
  //! given by \a flow_model for times between \a initial_time and \a final_time.
  SetModelType reachability_step(const FlowModelType& flow_model, 
                                 const SetModelType& initial_set_model, 
                                 const TimeType& initial_time, 
                                 const TimeType& final_time) const;
  
  //! \brief Computes the points reached by evolution of the \a initial_set_model under the flow
  //! given by \a flow_model for times between \a initial_time and \a final_time_model.
  SetModelType reachability_step(const FlowModelType& flow_model, 
                                 const SetModelType& initial_set_model, 
                                 const TimeType& initial_time, 
                                 const TimeModelType& final_time_model) const;
  
  //! \brief Computes the points reached by evolution of the \a initial_set_model under the flow
  //! given by \a flow_model for times between \a initial_time_model and \a final_time_model.
  SetModelType reachability_step(const FlowModelType& flow_model, 
                                 const SetModelType& initial_set_model, 
                                 const TimeModelType& initial_time_model, 
                                 const TimeModelType& final_time_model) const;
  
  //! \brief Gives the extended time model for the reachability step between the
  //! \a initial_time_model and the \a final_time_model. The new time is given by
  //! \f$\tau'(e,s) = (1-s)\tau_0(e)+s\tau_1(e)\f$.
  TimeModelType reachability_time(const TimeModelType& initial_time_model, 
                                  const TimeModelType& final_time_model) const;
  

  //! \brief A model for the map \a f over the domain \a d.
  MapModelType map_model(const FunctionType& f, BoxType& d) const;
  //! \brief A model for the flow determined by the vector field \a vf over the initial domain \a d,
  //! valid for times up to \a h, assuming that the state remains in the bounding box \a b.
  FlowModelType flow_model(const FunctionType& vf, const BoxType& d, 
                           const TimeType& h, const BoxType& b) const;
  //! \brief A model for the real-valued function \a g over the domain \a d.
  GuardModelType guard_model(const FunctionType& g, const BoxType& d) const;

  //! \brief Computed a pair \f$(h,B)\f$ such that the flow of the vector_field \a vf starting in
  //! domain \a d remains in \a B for times up to \a h. The maximum allowable \a h and maximum
  //! allowable diameter of \a B are given.
  pair<TimeType,BoxType> 
  flow_bounds(const FunctionType& vf, 
              const BoxType& d, 
              const RealType& maximum_step_size, 
              const RealType& maximum_bound_diameter) const;


  //@{ \name Conversion operations
  //! \brief Converts a set to a model (Unstable)
  SetModelType model(const SetType& s) const;
  //! \brief Converts a model to a set (Unstable)
  SetType set(const SetModelType& s) const;
  //@}

  //@{ \name Set-based operations
  //! \brief Tests if the set described by the model \a s is disjoint from the box \a box.
  tribool disjoint(const SetModelType& s, const BoxType& bx) const;
  //! \brief An upper bound for the radius of the set \a s.
  RealType radius(const SetModelType& s) const;
  //! \brief A box containing the set \a s.
  BoxType bounding_box(const SetModelType& s) const;
  //! \brief A list of sets obtained by subdividing the set \a s into at least two smaller pieces.
  array<SetModelType> subdivide(const SetModelType& s) const;
  //! \brief An over-approximation to the set \a s with a simplified description.
  SetModelType reduce(const SetModelType& s) const;
  //@}

};

}


#endif /* ARIADNE_DYNAMICAL_TOOLBOX_H */