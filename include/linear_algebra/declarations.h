/***************************************************************************
 *            linear_algebra/declarations.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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
 
/*! \file linear_algebra/declarations.h
 *  \brief Forward declarations of classes in the Linear Algebra module.
 */

#ifndef ARIADNE_LINEAR_ALGEBRA_DECLARATIONS_H
#define ARIADNE_LINEAR_ALGEBRA_DECLARATIONS_H


namespace Ariadne { 
  namespace LinearAlgebra {
    class MultiIndex;
      
    template<class R> class Vector;
    template<class R> class VectorSlice;
    template<class R> class Matrix;
    template<class R> class MatrixSlice;
    template<class R> class Tensor;
    template<class R> class LinearProgram;
  }
}

#endif /* ARIADNE_LINEAR_ALGEBRA_DECLARATIONS_H */