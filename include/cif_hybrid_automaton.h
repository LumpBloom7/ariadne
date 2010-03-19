/***************************************************************************
 *            cif_hybrid_automaton.h
 *
 *  Copyright  2009  Pieter Collins
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

/*! \file cif_hybrid_automaton.h
 *  \brief Hybrid system classes for interface with the Compositional Interchange Format.
 */

#ifndef ARIADNE_CIF_HYBRID_AUTOMATON_H
#define ARIADNE_CIF_HYBRID_AUTOMATON_H

#include "hybrid_automaton.h"


namespace Ariadne {



/*! \brief An interface to the Ariadne hybrid automaton class conforming to the Compositional Interchange Format (CIF)
 *  concept of an Atomic Interchange Automaton.
 *
 *  Currently, Ariadne does not support invariants, so only time-can-progress (tcp) predicates can be used.
 * \sa \link Ariadne::AtomicHybridAutomaton \endlink

 */
class CIFAtomicInterchangeAutomaton
    : public AtomicHybridAutomaton
{
  public:
    //@{
    //! \name Constructors and destructors

    //! \brief Construct an empty automaton with the given name
    CIFAtomicInterchangeAutomaton(const String& name);

    //@}

    //@{
    //! \name Methods for building the automaton.

    //! \brief Adds a discrete mode to the automaton, with given algebraic and differential equations and progress predicates.
    //!
    //!    \param location is the mode's location label.
    //!    \param algebraic_equations are the algebraic equations \f$y=f(x_1,\ldots,x_k)\f$ valid in the mode.
    //!    \param differential_equations are the differential equations \f$\dot{y}=f(x_1,\ldots,x_k)\f$ valid in the mode.
    //!    \param progress_predicates are the time-can-progress (tcp) predicates e.g. \f$g(x_1,\ldots,x_n)\leq c\f$ valid in the mode.
    //!
    //! Currently no invariant predicates are supported.
    void new_mode(AtomicDiscreteLocation location,
                  const List<RealAlgebraicAssignment>& algebraic_equations,
                  const List<RealDifferentialAssignment>& differential_equations,
                  const List<ContinuousPredicate>& progress_predicates);


    //! \brief Adds a discrete transition to the automaton using the discrete states to specify the source and target modes.
    //!
    //!    \param source is the transition's source location.
    //!    \param event is the transition's event.
    //!    \param target is the transition's target location.
    //!    \param reset is the transition's reset.
    //!    \param guard is the transition's activation region.
    //!    \param urgency is a flag giving whether the transition is urgent.
    void new_transition(AtomicDiscreteLocation source,
                        DiscreteEvent event,
                        AtomicDiscreteLocation target,
                        const List<RealUpdateAssignment>& reset,
                        const ContinuousPredicate& guard,
                        Urgency urgency=permissive);
};



inline
CIFAtomicInterchangeAutomaton::CIFAtomicInterchangeAutomaton(const String& name)
    : AtomicHybridAutomaton(name)
{
}

inline
void
CIFAtomicInterchangeAutomaton::new_mode(
    AtomicDiscreteLocation location,
    const List<RealAlgebraicAssignment>& algebraic_equations,
    const List<RealDifferentialAssignment>& differential_equations,
    const List<ContinuousPredicate>& progress_predicates)
{
    this->AtomicHybridAutomaton::new_mode(location,algebraic_equations,differential_equations);
    for(uint i=0; i!=progress_predicates.size(); ++i) {
        String name=String(location.name())+".tcp"+to_string(i);
        this->AtomicHybridAutomaton::new_invariant(location,DiscreteEvent(name),progress_predicates[i]);
    }
}

inline
void
CIFAtomicInterchangeAutomaton::new_transition(
    AtomicDiscreteLocation source,
    DiscreteEvent event,
    AtomicDiscreteLocation target,
    const List<RealUpdateAssignment>& reset,
    const ContinuousPredicate& guard,
    Urgency urgency)

{
    ARIADNE_ASSERT_MSG(urgency==permissive,"Only nonurgent transitions supported");
    this->AtomicHybridAutomaton::new_transition(source,event,guard,target,reset);
}




} // namespace Ariadne

#endif // ARIADNE_CIF_HYBRID_AUTOMATON_H
