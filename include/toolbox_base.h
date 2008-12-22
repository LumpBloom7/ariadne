/***************************************************************************
 *            toolbox_base.h
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
 
/*! \file toolbox_base.h
 *  \brief Base class for dynamical calculus routines containing implementations of routings which can be naturally composed from other routines. 
 */

#ifndef ARIADNE_TOOLBOX_BASE_H
#define ARIADNE_TOOLBOX_BASE_H

#include "tribool.h"
#include "logging.h"
#include "toolbox_interface.h"

/* \brief Top-level namespace. */
namespace Ariadne {
 
using std::pair;


template<class T> class array;

class Interval;
class FunctionInterface;
template<class X> class Vector;
class Box;



/*! \brief Tools for analysing dynamical systems based on function models. */
template<class Mdl> 
class ToolboxBase
    : public ToolboxInterface<Mdl>
    , public Loggable
{
    typedef Float R;
    typedef Float A;
    typedef Interval I;
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
    typedef Mdl PredicateModelType;
    typedef Float TimeType;

    typedef Vector<Interval> BoxType;
    typedef Interval IntervalType;
    typedef FunctionInterface FunctionType;
    typedef SetModelType EnclosureType;
  protected:
    tribool _tribool(const IntervalType& ivl) const { 
        if(ivl.lower()>0) { return true; } else if(ivl.upper()<0) { return false; } else { return indeterminate; } }
  public:
    //@{ \name Dynamical operations

    //! \brief Computes the image of the set defined by \a set_model under the approximation of the map 
    //! given by \a map_model.
    virtual SetModelType reset_step(const MapModelType& map_model, 
                                    const SetModelType& set_model) const = 0;
  
  
    //! \brief Test if a set satisfied the constraint given by the guard model. Returns \a true is all 
    //! points in the set satisfy the constraint, \a false if all points do not satisfy the constraint, 
    //! and indeterminate otherwise.
    virtual tribool active(const PredicateModelType& guard_model, 
                           const SetModelType& _set_model) const = 0;

    //! \brief Computes an over-approximation to the time interval for which the \a initial_set_model 
    //! touch the set specified by the \a guard model under the \a flow_model. The \a minimum and \a maximum_time 
    //! gives the minimum and maximum time for which the evolution is valid.
    virtual Interval touching_time_interval(const PredicateModelType& guard_model, 
                                            const FlowModelType& flow_model, 
                                            const SetModelType& initial_set_model) const = 0;

    //! \brief Computes the time at which points in the \a initial_set_model cross the zero-set of the
    //! the \a guard_model under evolution of the \a flow_model, for times between the \a minimum_time and \a maximum_time.
    //! The crossing must be (differentiably) transverse.
    virtual TimeModelType crossing_time(const PredicateModelType& guard_model,
                                        const FlowModelType& flow_model, 
                                        const SetModelType& initial_set_model) const = 0;

    //! \brief Computes the points reached by evolution of the \a initial_set_model under the flow
    //! given by \a flow_model. The \a integration_time_model \f$\tau(e)\f$ gives the time the point 
    //! starting at \f$x(e)\f$ should be flowed.
    virtual SetModelType integration_step(const FlowModelType& flow_model, 
                                          const SetModelType& initial_set_model, 
                                          const TimeModelType& integration_time_model) const = 0;
  
    //! \brief Computes the points reached by evolution of the \a initial_set_model under the flow
    //! given by \a flow_model for times between \a initial_time_model and \a final_time_model.
    virtual SetModelType reachability_step(const FlowModelType& flow_model, 
                                           const SetModelType& initial_set_model, 
                                           const TimeModelType& initial_time_model, 
                                           const TimeModelType& final_time_model) const = 0;
  
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
  
    //! \brief Computed a pair \f$(h,B)\f$ such that the flow of the vector_field \a vf starting in
    //! domain \a d remains in \a B for times up to \a h. The maximum allowable \a h and maximum
    //! allowable diameter of \a B are given.
    virtual pair<TimeType,BoxType> 
    flow_bounds(const FunctionType& vf, 
                const BoxType& d, 
                const RealType& maximum_step_size, 
                const RealType& maximum_bound_diameter) const = 0;

    //@}
   
    //@{ \name Constructing models for functions
    //! \brief A model for the map \a f over the domain \a d.
    virtual MapModelType map_model(const FunctionType& f, const BoxType& d) const = 0;

    //! \brief A model for the flow determined by the vector field \a vf over the initial domain \a d,
    //! valid for times up to \a h, assuming that the state remains in the bounding box \a b.
    virtual FlowModelType flow_model(const FunctionType& vf, const BoxType& d, 
                                     const TimeType& h, const BoxType& b) const = 0;

    //! \brief A model for the real-valued function \a g over the domain \a d.
    virtual PredicateModelType predicate_model(const FunctionType& g, const BoxType& d) const = 0;

    //! \brief A model for the constant time \a t over the box \a d.
    virtual TimeModelType time_model(const Float& t, const BoxType& d) const = 0;
    //@}

    //@{ \name Set-based operations
    //! \brief Compute a model for the given box \a bx.
    virtual SetModelType set_model(const BoxType& bx) const = 0;
    //! \brief Compute an enclosure for the set model \a s.
    virtual EnclosureType enclosure(const SetModelType& s) const = 0;
    //! \brief Tests if the set described by the model \a s is disjoint from the box \a box.
    virtual tribool disjoint(const SetModelType& s, const BoxType& bx) const = 0;
    //! \brief A box containing the set \a s.
    virtual BoxType bounding_box(const SetModelType& s) const = 0;
    //! \brief A list of sets obtained by subdividing the set \a s into at least two smaller pieces.
    virtual array<SetModelType> subdivide(const SetModelType& s) const = 0;
    //! \brief An over-approximation to the set \a s with a simplified description.
    virtual SetModelType simplify(const SetModelType& s) const = 0;
    //@}

  public:
    //! \brief Test if a box satisfies the constraint given by the guard. Returns \a true is all points
    //! in the box satisfy the constraint, \a false if all points do not satisfy the constraint, and 
    //! indeterminate otherwise.
    tribool active(const FunctionType& guard,  const BoxType& box) const { 
        return this->_tribool(guard.evaluate(box)[0]); }

    //! \brief Test if a set satisfied the constraint given by the guard. Returns \a true is all points
    //! in the set satisfy the constraint, \a false if all points do not satisfy the constraint, and 
    //! indeterminate otherwise.
    tribool active(const FunctionType& guard,  const SetModelType& set_model) const { 
        return this->active(this->predicate_model(guard,set_model.range()),set_model); }

 
    //! \brief Computes the image of the set defined by \a set_model under the \a map.
    SetModelType reset_step(const FunctionType& map, 
                            const SetModelType& set_model) const {
        return this->reset_step(this->map_model(map,set_model.range()),set_model); }
  
    //! \brief Computes the points reached by evolution of the \a initial_set_model under the flow
    //! given by \a flow_model. The \a integration_time gives the time all points should be flowed.
    SetModelType integration_step(const FlowModelType& flow_model, 
                                  const SetModelType& initial_set_model, 
                                  const TimeType& integration_time) const {
        return this->integration_step(flow_model,initial_set_model,this->time_model(integration_time,initial_set_model.domain())); }
  
    //! \brief Computes the points reached by evolution of the \a initial_set_model under the flow
    //! given by \a flow_model for times between \a initial_time and \a final_time.
    SetModelType reachability_step(const FlowModelType& flow_model, 
                                   const SetModelType& initial_set_model, 
                                   const TimeType& initial_time, 
                                   const TimeType& final_time) const {
        return this->reachability_step(flow_model,initial_set_model,this->time_model(initial_time,initial_set_model.domain()),this->time_model(final_time,initial_set_model.domain())); };
  
    //! \brief Computes the points reached by evolution of the \a initial_set_model under the flow
    //! given by \a flow_model for times between \a initial_time and \a final_time_model.
    SetModelType reachability_step(const FlowModelType& flow_model, 
                                   const SetModelType& initial_set_model, 
                                   const TimeType& initial_time, 
                                   const TimeModelType& final_time_model) const {
        return this->reachability_step(flow_model,initial_set_model,this->time_model(initial_time,initial_set_model.domain()),final_time_model); }
  
    //! \brief Computes the points reached by evolution of the \a initial_set_model under the flow
    //! given by \a flow_model for times between \a initial_time_model and \a final_time.
    SetModelType reachability_step(const FlowModelType& flow_model, 
                                   const SetModelType& initial_set_model, 
                                   const TimeModelType& initial_time_model, 
                                   const TimeType& final_time) const {
        return this->reachability_step(flow_model,initial_set_model,initial_time_model,this->time_model(final_time,initial_set_model.domain())); }
  
    //! \brief Gives the extended time model for the reachability step between the
    //! \a initial_time_model and the \a final_time_model. The new time is given by
    //! \f$\tau'(e,s) = (1-s)\tau_0(e)+s\tau_1(e)\f$.
    TimeModelType reachability_time(const TimeModelType& initial_time_model, 
                                    const TimeModelType& final_time_model) const;
  
    //! \brief Gives the extended time model for the reachability step between the
    //! \a initial_time and the \a final_time_model. The new time is given by
    //! \f$\tau'(e,s) = (1-s)\tau_0+s\tau_1(e)\f$.
    TimeModelType reachability_time(const TimeType& initial_time, 
                                    const TimeModelType& final_time_model) const {
        return this->reachability_time(this->time_model(initial_time,final_time_model.domain()),final_time_model); }
  
    //! \brief Gives the extended time model for the reachability step between the
    //! \a initial_time_model and the \a final_time. The new time is given by
    //! \f$\tau'(e,s) = (1-s)\tau_0(e)+s\tau_1\f$.
    TimeModelType reachability_time(const TimeModelType& initial_time_model, 
                                    const TimeType& final_time) const {
        return this->reachability_time(initial_time_model,this->time_model(final_time,initial_time_model.domain())); };
  
};


template<class Mdl>
typename ToolboxBase<Mdl>::TimeModelType
ToolboxBase<Mdl>::
reachability_time(const TimeModelType& initial_time_model, 
                  const TimeModelType& final_time_model) const
{
    ARIADNE_ASSERT(initial_time_model.argument_size()==final_time_model.argument_size());
    uint ng=initial_time_model.argument_size();

    ModelType expanded_initial_time_model=embed(initial_time_model,Vector<I>(ng+1,I(-1,1)),Vector<R>(ng+1,R(0)),0u);
    ModelType expanded_final_time_model=embed(final_time_model,Vector<I>(ng+1,I(-1,1)),Vector<R>(ng+1,R(0)),0u);
  
    ModelType time_interval_model=ModelType::affine(I(-1,1),R(0),A(0.5),A(0.5),1u,0u);
    ModelType expanded_time_interval_model=embed(time_interval_model,Vector<I>(ng+1,I(-1,1)),Vector<R>(ng+1,R(0)),ng);
    ModelType expanded_reach_time_model=expanded_initial_time_model+expanded_time_interval_model*(expanded_final_time_model-expanded_initial_time_model);

    return expanded_reach_time_model;
}



}




#endif /* ARIADNE_DYNAMICAL_TOOLBOX_H */
