/***************************************************************************
 *            hybrid_evolution_methods.h
 *
 *  Copyright  2004-7  Pieter Collins
 *  Pieter.Collins@cwi.nl
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

/*! 

\file hybrid_evolution_methods.h
\brief Documentation on methods for evolution of hybrid systems



\page hybrid_evolution Hybrid Evolution Methods

The evolution of a hybrid system takes place in a number of steps. 
Each step consists of a continuous evolution, followed by at most one discrete transition. 
Since the time of a discrete transition may depend on the initial state, the evolution time needs to be stored as well as the point. Both the set and evolution times are stored as a polynomial "model".


\rationale

An alternative approach is to take evolution steps of a fixed time interval. However, this means that a point may undergo several discrete events in an evolution step, and the dispatching of this evolution may be hard to compute. 


\section hybrid_timed_set Timed hybrid sets


Due to the need to keep track of transition times, a hybrid evolution is defined on timed hybrid sets. 
A timed hybrid set consists of a model for the state and for the evolution time. The evolution time is either a rational constant or a function model (e.g. an affine model) defined on the same variables as the set model: \f[ X = c + Ge; \quad t=s+re. \f]
Where possible, an evolution will end at an exact time. 

\section evolution_traces Evolution traces

An evolution dispatcher is an algorithm which computes how long the integration needs to procees, and which discrete events need to be processed. 
The maximum allowable integration time is computed by evaluating all the constraint conditions and testing for the smallest integration time.
If the flow is transverse to the constraint set, then a differentiable model for the integration time can be constructed.

First, the starting set is evalutated to test if any events are (partially) activated or any invariants (partially) not satisfied.
 - If a guard is totally activated, then the event occurs immediately.
 - If an overflowing guard is partially activated, then the evolution is computed with both the event occurring immediately, and a flow taking place. The initial set may first be subdivided.
 - If a guard is partially activated but the set is repelling, then both an immediate event and a flow occur. The initial set may first be subdivided.

If the constraint with the smallest integration time is indeterminate, or depends on the continuous state, then either an upper bound is taken, or the set is subdivided.

The integration ends at a fixed time
 - If the maximum integration time is reached, or
 - If no events are activated.

The box may be subdivided if
 - The crossing is transverse, but takes longer than half the integration step.
The integration time is reduced if
 - The crossing is transverse, is started after half the integration step, but is not completed within the integration step.
 
If one of the smallest integration times corresponds to a guard, then the guard is activated and the flow to the guard set is computed.



\subsection upper_evolution_trace Upper evolution traces

 -# Compute a bound \f$B\f$ for the flow \f$\Phi(X,[0,h])\f$.
 -# For each constraint \f$g_e\f$:
     - Compute an approximation to \f$g(B)\f$ and determine whether the constraint is satisfied, unsatisfied or crossed.
     - For each constraint which is crossed, estimate the switching time \f$s_e(x)\f$:
        - If the crossing is transverse, give a first-order approximation to \f$s_e(x)\f$.
        - If the crossing is not transverse, give a constant lower bound for the crossing time.
 -# If \f$s_e\f$ becomes negative for some constraint \f$e\f$:
     - If the radius of \f$X\f$ exceeds \c maximum_splitting_set_radius, subdivide \f$X\f$
     - Otherwise, perform both an \f$e\f$ and \f$t\f$ step.
 -# Compute the maximum flow time \f$\tau(x)\f$ and discard all events whose time exceeds the maximum flow time.
 -# If more than one blocking event is active, and the radius of \f$X\f$ exceeds the maximum_splitting_set_radius, subdivide \f$X\f$.
 -# If only one blocking event is active, and the crossing time interval exceeds the maximum_crossing_time, subdivide \f$X\f$.

\subsection lower_evolution_trace Lower evolution traces

In a lower evolution trace, if a box needs to be split, then the evolution is terminated, unless
 - The split occurs due to an unforced event, in which case both integration and an event occur
 - The split occurs due to mapping to a guard set which is overflowing, in which case both events occur, and the C<sup>0</sup> union of the resulting boxes is taken.
The evolution is also terminated if the computation cannot determine whether an activation or an invariant is first crossed, or which of two guards is first crossed.

\section switching_time Switching time

Let \f$s(x)\f$ be the time needed to flow from \f$x\f$ to the guard set \f$g(x)=0\f$.
The switching time \f$s(x)\f$ satisfies the equation \f$g\bigr(\Phi_1(x,s(x))\bigl) = 0\f$.
Then \f$s(x)\f$ satisfies
\f[ - \nabla s(x) = \frac{\nabla g(\Phi_1(x,s(x)))\cdot D\Phi_1(x,s(x))}{\nabla g(\Phi_1(x,s(x))\cdot f_1(\Phi_1(x,s(x)))} 
                  = \frac{\nabla g(y)\cdot D\Phi_1(x,s(x))}{\nabla g(y)\cdot f_1(y)} 
\f]

\section forced_transitions Forced transitions

Let \f$s(x)\f$ be the time needed to flow from \f$x\f$ to the guard set \f$g(x)=0\f$.
Then the transition is given by 
\f[ \Psi(x,t) = \Phi_2(r(\Phi_1(x,s(x))),t-s(x)) . \f]
and the Jacobian derivative is
\f[ D\Psi(x,t) \in D\Phi_2(B_2)\,Dr(B_1)\,D\Phi_1(B_1) \, + \, \bigl( D\Phi_2(B_2) \, Dr(B_1) \, f_1(B_1) - f_2(B_2) \bigr) \nabla s(x) . \f]
Suppose \f$s(c)\f$ is known. Then
\f[ \Psi(x,t) = \Psi(c,t) + D\Phi_2(B_2) Dr(B_1) D\Phi_1(B_1) + (f_1(B_1) - f_2(B_2)) \nabla\tau(x) r(\Phi_1(x,\tau(x))),t-\tau(x)) . \f]

The <em>saltation map</em> \f$\Psi(x,0)\f$ satisfies
\f[\begin{aligned}
    \Psi(x,0) &:= \Phi_2\bigl(r(\Phi_1(x,s(x))),-s(x)\bigr) \\
            &=  \bigr[ \Phi_2\bigl(r(\Phi_1(x,s(x))),-s(x)\bigr) - \Phi_2\bigl(r(\Phi_1(x,s(x))),0\bigr) \bigr] + \Phi_2\bigl(r(\Phi_1(x,s(x))),0\bigr) \\
            &= -\dot{\Phi}_2\bigl(r(\Phi_1(x,s(x))),\sigma\bigr)\,s(x) + r(\Phi_1(x,s(x))) \\
            &= -f_2\bigl(\Phi_2\bigl(r(\Phi_1(x,s(x))),\sigma\bigr)\bigr)\,s(x) + \bigr[ r(\Phi_1(x,s(x))) - r(\Phi_1(x,0)) \bigr] + r(\Phi_1(x,0)) \\
            &= -f_2(\zeta)\, s(x) + Dr(\Phi_1(x,\tau))\,f_1(\Phi_1(x,\tau))\,s(x) + r(x) \\
            &= -f_2(\zeta)\,s(x) + Dr(\eta)\,f(\eta)\,s(x) + r(c) + Dr(\xi)\,(x-c) \\[\jot]

            &\in r(c) + Dr(\xi)\cdot(x-c) + \bigl(Dr(B_1)\,f(B_1)-f_2(B_2)\bigr)\cdot s(x)
\end{aligned}\f] 

The evolution \f$\Psi(x,t)\f$ satisfies
\f[\begin{aligned}
    \Psi(x,t) &:= \Phi_2\bigl(r(\Phi_1(x,s(x))),t-s(x)\bigr) \\
            &=  \bigr[ \Phi_2\bigl(r(\Phi_1(x,s(x))),t-s(x)\bigr) - \Phi_2\bigl(r(\Phi_1(x,s(x))),t-s(c)\bigr) \bigr] + \Phi_2\bigl(r(\Phi_1(x,s(x))),t-s(c)\bigr) \\
            &= \dot{\Phi}_2\bigl(r(\Phi_1(x,s(x))),\sigma\bigr)\,(s(c)-s(x)) + \Phi_2\bigl(r(\Phi_1(x,s(x))),t-s(c)\bigr) \\
            &= -f_2\bigl(\Phi_2\bigl(r(\Phi_1(x,s(x))),\sigma\bigr)\bigr)\,(s(c)-s(x)) + \bigr[ \Phi_2\bigl(r(\Phi_1(x,s(x))),t-s(c)\bigr) - \Phi_2\bigl(r(\Phi_1(x,s(c))),t-s(c)\bigr) \bigr] + \Phi_2\bigl(r(\Phi_1(x,s(c))),t-s(c)\bigr) \\
            &= -f_2(\zeta)\, s(x) + D\Phi_2\bigl(r(\Phi_1(x,\tau)),t-s(c)\bigr) Dr(\Phi_1(x,\tau))\,f_1(\Phi_1(x,\tau))\,(s(x)-s(c)) + \Phi_2\bigl(r(\Phi_1(x,s(c))),t-s(c)\bigr)  \\
\end{aligned}\f] 
If \f$s\f$ is an approximation to \f$s(c)\f$, then the evolution \f$\Psi(x,t)\f$ can be computed by
\f[\begin{aligned}
    \Psi(x,t) &= \Phi_2\bigl(\Psi(\Phi_1(x,s),0),t-s\bigr) = \Phi_2^{t-s} \circ \Psi^0 \circ \Phi_1^s(x,t) 
\end{aligned}\f] 



\subsection nonsmooth_forced_transitions Nonsmooth forced transitions

If the switching time \f$s(x)\f$ is non-smooth or discontinuous, we need to rely on zero-order methods for computing the evolution.
If the switching time is bounded by \f$[-h,+h]\f$, then the saltation map is given by
\f[\begin{aligned}
    \Psi(x,0) &:= \Phi_2\bigl(r(\Phi_1(x,s(x))),-s(x)\bigr) \\[\jot]
              &\in  r(X+[-h,h]f_1(B_1))+[-h,h]f_2(B_2)
\end{aligned}\f] 


\section unforced_transitions Unforced transitions

Consider a transition which can occur at any time in the interval \f$[-h,+h\f$]. Suppose \f$X\f$ is a bound for \f$x\f$, \f$B_1\f$ is a bound for \f$\Phi_1(X,[-h,h])\f$ and \f$B_2\f$ is a bound for \f$\Phi_2\bigl(r(\Phi_1(X,t)),-t\bigr)\f$ for \f$t\in[-h,+h]\f$. We have
\f[ \begin{aligned} \Psi(x,0;t) &:= \Phi_2(r(\Phi_1(x,t)),-t) \\
                                &=  \bigl[ \strut \Phi_2(r(\Phi_1(x,t)),-t) - \Phi_2(r(\Phi_1(x,t),0)) \bigr] + \Phi_2(r(\Phi_1(x,t),0)) \\
                                &=  -t\,\dot{\Phi}_2\bigl(r(\Phi_1(x,t)),\tau\bigr) + r(\Phi_1(x,t)) \\
                                &=  -t\,f_2\left(\Phi_2\bigl(r(\Phi_1(x,t)),\tau\bigr)\right) + \bigl[ r(\Phi_1(x,t)) - r(\Phi_1(x,0)) \bigr] + r(\Phi_1(x,0)) \\
                                &= -t\,f_2(\zeta) + t\,Dr(\Phi_1(x,\tau))\,f_1(\Phi_1(x,\tau)) + r(x) \\
                                &= \bigl(Dr(\eta)\,f_1(\eta)-f_2(\zeta)\bigr)\cdot t + Dr(\xi)\cdot(x-c) + r(c) \\[\jot]
                                &\in r(c) + Dr(X)\cdot(x-c) + \bigl(Dr(B_1)\,f_1(B_1) - f_2(B_2) \bigr) \cdot t 
\end{aligned} \f]
If we additionally wish to flow forward a time \f$h\f$, then 
*/