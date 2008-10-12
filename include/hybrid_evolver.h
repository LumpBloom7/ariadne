/***************************************************************************
 *            hybrid_evolver.h
 *
 *  Copyright  2007-8  Alberto Casagrande, Pieter Collins
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
 
/*! \file hybrid_evolver.h
 *  \brief Evolver for hybrid systems.
 */

#ifndef ARIADNE_HYBRID_EVOLVER_H
#define ARIADNE_HYBRID_EVOLVER_H

#include <string>
#include <vector>
#include <list>
#include <iostream>

#include <boost/smart_ptr.hpp>

#include "tuple.h"

#include "hybrid_automaton.h"
#include "evolver_interface.h"
#include "evolver_base.h"
#include "evolution_parameters.h"

namespace Ariadne {  
  
template<class Sys, class BS> class Evolver;

class ApproximateTaylorModel;
class HybridAutomaton;

class EvolutionParameters;
template<class MDL> class DynamicalToolbox;

class EvolutionProfiler;

typedef std::pair<int,Float> HybridTime;


typedef ApproximateTaylorModel DefaultModelType;

/*! \ingroup Evolve 
 *  \brief A class for computing the evolution of a hybrid system. 
 *
 * The actual evolution steps are performed by the HybridEvolver class.
 */
class HybridEvolver
  : public EvolverBase< HybridAutomaton, DefaultModelType >
{
  typedef Ariadne::DefaultModelType ModelType;
 public:
  typedef Float RealType;
  typedef HybridAutomaton SystemType;
  typedef ModelType ContinuousEnclosureType;
  typedef pair<DiscreteState,ContinuousEnclosureType> EnclosureType;
  typedef ListSet<EnclosureType> EnclosureListType;
  typedef HybridTime TimeType;
  typedef Float ContinuousTimeType;
 public:
    
  //! \brief Default constructor.
  HybridEvolver();
  
  //! \brief Construct from parameters using a default integrator.
  HybridEvolver(const EvolutionParameters& parameters);
  

  //@{
  //! \name Parameters controlling the evolution.
  //! \brief A reference to the parameters controlling the evolution.
  EvolutionParameters& parameters() { return *this->_parameters; }
  const EvolutionParameters& parameters() const { return *this->_parameters; }

  //@}
  

  //@{
  //! \name Evolution using abstract sets.
  //! \brief Compute an approximation to the evolution set using upper semantics. 
  EnclosureListType evolve(const SystemType& system, const EnclosureType& initial_set, const TimeType& time) const {
    EnclosureListType final; EnclosureListType reachable; EnclosureListType intermediate; 
    this->_evolution(final,reachable,intermediate,system,initial_set,time,upper_semantics,false); 
        return final; }

  //! \brief Compute an approximation to the evolution set under upper semantics. 
  EnclosureListType reach(const SystemType& system, const EnclosureType& initial_set, const TimeType& time) const {
        EnclosureListType final; EnclosureListType reachable; EnclosureListType intermediate; 
        this->_evolution(final,reachable,intermediate,system,initial_set,time,upper_semantics,true); 
        return intermediate; }

 protected:
  virtual void _evolution(EnclosureListType& final, EnclosureListType& reachable, EnclosureListType& intermediate, 
                          const SystemType& system, const EnclosureType& initial, const TimeType& time, 
                          Semantics semantics, bool reach) const;

 private:
  boost::shared_ptr< EvolutionParameters > _parameters;
  boost::shared_ptr< DynamicalToolbox<ModelType> > _toolbox;
  //boost::shared_ptr< EvolutionProfiler >  _profiler;
};


  
} // namespace Ariadne

#endif // ARIADNE_HYBRID_EVOLVER_H
