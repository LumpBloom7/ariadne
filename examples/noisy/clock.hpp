/***************************************************************************
 *            clock.hpp
 *
 *  Copyright  2008-18 Luca Geretti
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

#include "ariadne.hpp"

using namespace Ariadne;


Tuple<String,DottedRealAssignments,RealVariablesBox,RealVariablesBox,Real,double> CK()
{
    RealVariable x("x"), y("y"), u1("u1"), u2("u2");
    DottedRealAssignments dynamics={dot(x)=u1,dot(y)=u2};
    RealVariablesBox inputs={0.95_dec<=u1<=1.05_dec,0.95_dec<=u2<=1.05_dec};

    Real e=1/128_q;
    RealVariablesBox initial={{x==0},{y==0}};

    auto evolution_time=5;
    double step=1.0/256;

    return make_tuple("CK",dynamics,inputs,initial,evolution_time,step);
}
