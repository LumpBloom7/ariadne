/***************************************************************************
 *            python/float.h
 *
 *  Copyright  2005  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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

/*! \file python/float.h
 *  Typedef for Python floating-point class.
 */
 
#ifndef ARIADNE_PYTHON_FLOAT_H
#define ARIADNE_PYTHON_FLOAT_H

#include <config.h>

#include "numeric/float.template.h"
#include "numeric/float.inline.h"

#if PYTHON_FLOAT == Float64 

#include "numeric/float64.h" 
namespace Ariadne { 
  namespace Python { 
    typedef Numeric::Float64 FloatPy; 
    typedef Numeric::Interval<FloatPy> IntervalPy;
  } 
}

#elif PYTHON_FLOAT == FloatMP 
#include "numeric/floatmp.h" 
namespace Ariadne { 
  namespace Python { 
    typedef Numeric::FloatMP FloatPy; 
    typedef Numeric::Interval<FloatPy> IntervalPy;
  } 
}

#endif

#endif /* ARIADNE_PYTHON_FLOAT_H */