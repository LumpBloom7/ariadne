/***************************************************************************
 *            valve-permissive.hpp
 *
 *  Copyright  2017 Luca Geretti
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

#include <cstdarg>
#include "ariadne.hpp"

using namespace Ariadne;

HybridAutomaton getValve()
{

    // Declare some constants. Note that system parameters should be given as variables.
    RealConstant T("T",4.0_decimal);
    RealConstant hmin("hmin",5.5_decimal);
    RealConstant hmax("hmax",8.0_decimal);
    RealConstant delta("delta",0.05_decimal);

    // Declare the shared system variables
    RealVariable aperture("aperture");
    RealVariable height("height");

    // Declare the events we use
    DiscreteEvent start_opening("start_opening");
    DiscreteEvent start_closing("start_closing");
    DiscreteEvent finished_opening("finished_opening");
    DiscreteEvent finished_closing("finished_closing");
    DiscreteEvent must_start_opening("must_start_opening");
    DiscreteEvent must_start_closing("must_start_closing");

    HybridAutomaton valve_automaton("valve_automaton");
    StringVariable valve("valve");

    // Declare the values the valve can variable can have
    StringConstant open("open");
    StringConstant opening("opening");
    StringConstant closed("closed");
    StringConstant closing("closing");

    // Since aperture is a known constant when the valve is open or closed,
    // specify aperture by an algebraic equation.
    valve_automaton.new_mode(valve|open,{let(aperture)=+1.0_decimal});
    valve_automaton.new_mode(valve|closed,{let(aperture)=0.0_decimal});
    // Specify the differential equation for how the valve opens/closes.
    valve_automaton.new_mode(valve|opening,{dot(aperture)=+1/T});
    valve_automaton.new_mode(valve|closing,{dot(aperture)=-1/T});

    // Specify the invariants valid in each mode. Note that every invariant
    // must have an action label. This is used internally, for example, to
    // check non-blockingness of urgent actions.
    valve_automaton.new_invariant(valve|open,height<=hmax+delta,must_start_closing);
    valve_automaton.new_invariant(valve|closed,height>=hmin-delta,must_start_opening);

    valve_automaton.new_transition(valve|closed,start_opening,valve|opening,{next(aperture)=aperture},height<=hmin+delta,permissive);
    valve_automaton.new_transition(valve|opening,finished_opening,valve|open,aperture>=1,urgent);
    valve_automaton.new_transition(valve|open,start_closing,valve|closing,{next(aperture)=aperture},height>=hmax-delta,permissive);
    valve_automaton.new_transition(valve|closing,finished_closing,valve|closed,aperture<=0,urgent);

    return valve_automaton;
}
