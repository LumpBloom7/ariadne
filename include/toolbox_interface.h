/***************************************************************************
 *            toolbox_interface.h
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
 
/*! \file toolbox_interface.h
 *  \brief Interfaces for tools for working with dynamical systems. 
 */


#ifndef ARIADNE_TOOLBOX_INTERFACE_H
#define ARIADNE_TOOLBOX_INTERFACE_H

#include "tribool.h"

namespace Ariadne {
 
template<class T> class array;

class Interval;
class FunctionInterface;
template<class X> class Vector;



/*! \brief Tools for analysing dynamical systems based on function models. */
template<class Mdl> 
class ToolboxInterface
{
    typedef Float R;
    typedef Float A;
    typedef Interval I;
  private:
    ushort _spacial_order;
    ushort _temporal_order;
    ushort _order;
    ushort _smoothness;
  public:
    //!
    typedef Float RealType;
    //!
    typedef Mdl ModelType;
    typedef Mdl SetType;
    //!
    typedef Mdl SetModelType;
    typedef Mdl TimeModelType;
    typedef Mdl MapModelType;
    typedef Mdl FlowModelType;
    typedef Mdl GuardModelType;
    typedef Float TimeType;
    //typedef Box<RealType> BoxType;
    typedef Vector<Interval> BoxType;
    typedef Interval IntervalType;
    typedef FunctionInterface FunctionType;
  public:
    //! \brief Virtual destructor.
    virtual ~ToolboxInterface() { }

    //! \brief Compose a function with a set model.
    virtual SetModelType 
    compose(const FunctionInterface& function, 
            const SetModelType& set_model) const = 0;
  
    //! \brief Compose a function model with a set model.
    virtual SetModelType 
    compose(const FlowModelType& function_model, 
            const SetModelType& set_model) const = 0;
  
    //! \brief Test if a box satisfies the constraint given by the guard. Returns \a true is all points
    //! in the box satisfy the constraint, \a false if all points do not satisfy the constraint, and 
    //! indeterminate otherwise.
    virtual tribool 
    active(const FunctionType& guard,
           const BoxType& box) const = 0;

    //! \brief Test if a set satisfied the constraint given by the guard. Returns \a true is all points
    //! in the set satisfy the constraint, \a false if all points do not satisfy the constraint, and 
    //! indeterminate otherwise.
    virtual tribool 
    active(const FunctionType& guard,  
           const SetModelType& set_model) const = 0;

    //! \brief Test if a set satisfied the constraint given by the guard model. Returns \a true is all 
    //! points in the set satisfy the constraint, \a false if all points do not satisfy the constraint, 
    //! and indeterminate otherwise.
    virtual tribool 
    active(const GuardModelType& guard_model, 
           const SetModelType& _set_model) const = 0;

    //! \brief Computes an over-approximation to the time interval for which the \a initial_set_model 
    //! touch the set specified by the \a guard model under the \a flow_model. The \a minimum and \a maximum_time 
    //! gives the minimum and maximum time for which the evolution is valid.
    virtual Interval
    touching_time_interval(const GuardModelType& guard_model, 
                           const FlowModelType& flow_model, 
                           const SetModelType& initial_set_model) const = 0;  

    //! \brief Computes the time at which points in the \a initial_set_model cross the zero-set of the
    //! the \a guard_model under evolution of the \a flow_model, for times between the \a minimum_time and \a maximum_time.
    //! The crossing must be (differentiably) transverse.
    virtual TimeModelType 
    crossing_time(const FlowModelType& flow_model,
                  const GuardModelType& guard_model,
                  const SetModelType& initial_set_model) const = 0;
  
    //! \brief Computes the image of the set defined by \a set_model under the \a map.
    /*
      virtual SetModelType 
      reset_step(const FunctionType& map, 
      const SetModelType& set_model) const = 0;  
    */

    //! \brief Computes the image of the set defined by \a set_model under the approximation of the map 
    //! given by \a map_model.
    /*
      virtual SetModelType 
      reset_step(const MapModelType& map_model, 
      const SetModelType& set_model) const = 0;
    */

    //! \brief Computes the points reached by evolution of the \a initial_set_model under the flow
    //! given by \a flow_model. The \a integration_time gives the time all points should be flowed.
    virtual SetModelType 
    integration_step(const FlowModelType& flow_model, 
                     const SetModelType& initial_set_model, 
                     const TimeType& integration_time) const = 0;
  
    //! \brief Computes the points reached by evolution of the \a initial_set_model under the flow
    //! given by \a flow_model. The \a integration_time_model \f$\tau(e)\f$ gives the time the point 
    //! starting at \f$x(e)\f$ should be flowed.
    virtual SetModelType 
    integration_step(const FlowModelType& flow_model, 
                     const SetModelType& initial_set_model, 
                     const TimeModelType& integration_time_model) const = 0;
  
    //! \brief Computes the points reached by evolution of the \a initial_set_model under the flow
    //! given by \a flow_model for times given by \a reachability_time_model. 
    //! The \a reachability_time_model must have one more independent variable than the 
    //! \a initial_set_model.
    //! 
    //! \invariant <code>reachability_time_model.argument_size()==initial_set_model.argument_size()+1</code>
    virtual SetModelType 
    reachability_step(const FlowModelType& flow_model, 
                      const SetModelType& initial_set_model, 
                      const TimeModelType& reachability_time_model) const = 0;
  
    //! \brief Computes the points reached by evolution of the \a initial_set_model under the flow
    //! given by \a flow_model for times between \a initial_time and \a final_time.
    virtual SetModelType 
    reachability_step(const FlowModelType& flow_model, 
                      const SetModelType& initial_set_model, 
                      const TimeType& initial_time, 
                      const TimeType& final_time) const = 0;
  
    //! \brief Computes the points reached by evolution of the \a initial_set_model under the flow
    //! given by \a flow_model for times between \a initial_time and \a final_time.
    virtual SetModelType 
    reachability_step(const FlowModelType& flow_model, 
                      const SetModelType& initial_set_model, 
                      const TimeModelType& initial_time_model, 
                      const TimeType& final_time) const = 0;
  
    //! \brief Computes the points reached by evolution of the \a initial_set_model under the flow
    //! given by \a flow_model for times between \a initial_time and \a final_time_model.
    virtual SetModelType 
    reachability_step(const FlowModelType& flow_model, 
                      const SetModelType& initial_set_model, 
                      const TimeType& initial_time, 
                      const TimeModelType& final_time_model) const = 0;
  
    //! \brief Computes the points reached by evolution of the \a initial_set_model under the flow
    //! given by \a flow_model for times between \a initial_time_model and \a final_time_model.
    virtual SetModelType
    reachability_step(const FlowModelType& flow_model, 
                      const SetModelType& initial_set_model, 
                      const TimeModelType& initial_time_model, 
                      const TimeModelType& final_time_model) const = 0;
  
    //! \brief Gives the extended time model for the reachability step between the
    //! \a initial_time_model and the \a final_time_model. The new time is given by
    //! \f$\tau'(e,s) = (1-s)\tau_0(e)+s\tau_1(e)\f$.
    virtual TimeModelType
    reachability_time(const TimeModelType& initial_time_model, 
                      const TimeModelType& final_time_model) const = 0;
  
    //! \brief Gives the extended time model for the reachability step between the
    //! \a initial_time and the \a final_time_model. The new time is given by
    //! \f$\tau'(e,s) = (1-s)\tau_0+s\tau_1(e)\f$.
    virtual TimeModelType
    reachability_time(const TimeType& initial_time, 
                      const TimeModelType& final_time_model) const = 0;
  
    //! \brief Gives the extended time model for the reachability step between the
    //! \a initial_time_model and the \a final_time. The new time is given by
    //! \f$\tau'(e,s) = (1-s)\tau_0(e)+s\tau_1\f$.
    virtual TimeModelType
    reachability_time(const TimeModelType& initial_time_model, 
                      const TimeType& final_time) const = 0;
  
  
    //! \brief A model for the map \a f over the domain \a d.
    virtual MapModelType
    map_model(const FunctionType& f, 
              const BoxType& d) const = 0;

    //! \brief A model for the flow determined by the vector field \a vf over the initial domain \a d,
    //! valid for times up to \a h, assuming that the state remains in the bounding box \a b.
    virtual FlowModelType
    flow_model(const FunctionType& vf, 
               const BoxType& d, 
               const TimeType& h, 
               const BoxType& b) const = 0;
  
    //! \brief A model for the real-valued function \a g over the domain \a d.
    virtual GuardModelType 
    predicate_model(const FunctionType& g, 
                    const BoxType& d) const = 0;

    //! \brief A model for the constant time function \a t over the domain \a d.
    virtual TimeModelType 
    constant_time_model(const Float& t, 
                        const BoxType& d) const = 0;

  
    //! \brief Computed a pair \f$(h,B)\f$ such that the flow of the vector_field \a vf starting in
    //! domain \a d remains in \a B for times up to \a h. The maximum allowable \a h and maximum
    //! allowable diameter of \a B are given.
    virtual std::pair<TimeType,BoxType> 
    flow_bounds(const FunctionType& vf, 
                const BoxType& d, 
                const RealType& maximum_step_size, 
                const RealType& maximum_bound_diameter) const = 0;


    //@{ \name Conversion operations

    //! \brief Converts a set to a model (Unstable)
    virtual SetModelType model(const SetType& s) const = 0;

    //! \brief Converts a model to a set (Unstable)
    virtual SetType set(const SetModelType& s) const = 0;

    //@}

    //@{ \name Set-based operations
    //! \brief Tests if the set described by the model \a s is disjoint from the box \a box.
    virtual tribool disjoint(const SetModelType& s, const BoxType& bx) const = 0;
    //! \brief A box containing the set \a s.
    virtual BoxType bounding_box(const SetModelType& s) const = 0;
    //! \brief A list of sets obtained by subdividing the set \a s into at least two smaller pieces.
    virtual array<SetModelType> subdivide(const SetModelType& s) const = 0;
    //! \brief An over-approximation to the set \a s with a simplified description.
    virtual SetModelType simplify(const SetModelType& s) const = 0;
    //@}

};

} //  namespace Ariadne


#endif // ARIADNE_TOOLBOX_INTERFACE_H */
