/***************************************************************************
 *            function/identity_function.cc
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
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

#include "numeric/float.h"

#include "function/identity_function.h"
#include "function/identity_function.code.h"

namespace Ariadne {
  namespace Function {
    using namespace Numeric;
    
    template class IdentityFunction<Rational>;

#ifdef ENABLE_FLOAT64
    template class IdentityFunction<Float64>;
#endif
    
#ifdef ENABLE_FLOATMP
    template class IdentityFunction<FloatMP>;
#endif

  }
}