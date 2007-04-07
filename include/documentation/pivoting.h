/***************************************************************************
 *            pivoting.h
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

\file pivoting.h
\brief Documentation on pivoting and permutations



\page pivot Pivoting and Permutations

Re-ordering of matrix rows/columns and vector elements can be described by <em>permutations</em> and <em>pivoting</em>.
Essentially, a permutation gives the re-ordering directly, and is most useful if the re-ordered elements are to be used in-place or copied, whereas pivoting describes the re-ordering as a product of transpositions, and is most useful if the elements are to be re-ordered in place.

We write a permuation
\f[ \pi=\left(\begin{array}{cccccc}0&1&\cdots&n\!-\!1\\\pi_0&\pi_1&\cdots&\pi_{n-1}\end{array}\right); \quad \pi_i\neq\pi_j \text{ if } i\neq j \f]
and a pivot
\f[ p = (0\;p_0)\ (1\;p_1)\ \cdots\ (n\!-\!1,p_{n-1}); \quad p_i\geq i . \f]
where the transpositions are applied from left to right.

For example, the permutation
\f[ \pi=\left(\begin{array}{ccccc}0&1&2&3&4\\1&3&4&2&0\end{array}\right)  \f]
can be writing in pivot form
\f[ p=(0\;1)\ (1\;3)\ (2\;4)\ (3\;4)\ (4\;4) . \f]
The inverses are
\f[ \pi^{-1}=\left(\begin{array}{ccccc}0&1&2&3&4\\4&0&3&1&2\end{array}\right); \qquad 
     p^{-1}=(4\;4)\ (3\;4)\ (2\;4)\ (1\;3)\ (0\;1) \f]

Applying the permuation to the vector \f$v=(v_0,v_1,v_2,v_3,v_4)^T\f$ yields the vector \f$\pi v=(v_1,v_3,v_4,v_2,v_0)^T\f$.
Applying the pivots to the vector yield successively \f$(v_1,v_0,v_2,v_3,v_4)^T\f$, \f$(v_1,v_3,v_2,v_0,v_4)^T\f$, \f$(v_1,v_3,v_4,v_0,v_2)^T\f$ and \f$Pv=(v_1,v_3,v_4,v_2,v_0)^T\f$.

The pivot form is used by linear algebra packages such as LAPACK.
For the triangular factorization of a matrix \f$A\f$ with partial pivoting yields matrices \f$P,L,U\f$ such that \f$PA=LU\f$, where \f$P\f$ is a pivot matrix.
The matrix
\f[ A=\left(\begin{matrix}0&0&0&0&1\\1&x&x&x&x\\0&0&0&1&x\\0&1&x&x&x\\0&0&1&x&x\end{matrix}\right)\f]
yields a pivot matrix \f$P\f$ described by the pivot \f$p=(0\;1)\ (1\;3)\ (2\;4)\ (3\;4)\ (4\;4)\f$. In other words, during the factorisation, we first interchange row 0 and row 1, then row 1 (which corresponds to the old row 0) with row 3, and so on.

The permuation form is mostly used by linear programming packages to select basis variables.
<br>
Consider the \f$m\times n\f$ matrix \f$\mathrm{A}\f$ with \f$n>m\f$, and suppose \f$\mathrm{A}\f$ has full row rank.
Let \f$B=\{ j_0,j_1,\ldots,j_{m-1}\}\f$ be a set of \f$m\f$ values from \f$\{0,\ldots,n\!-\!1\}\f$, and define the matrix \f$\mathrm{B}=\mathrm{A}_B=\bigl( a_{i_0}\ a_{i_1}\ \cdots\ a_{i_{n-1}}\bigr)\f$.
Here, a (partial) permutation is used to select the pivot elements.
<br>
In the reduced simplex algorithm, we store an intermediate form of the problem as \f$x_B+\mathrm{A}x_N=b\f$
We now store both \f$B=(j_0,j_1,\ldots,j_{m-1})\f$ giving the columns of the basic variables, and \f$N=(j_{m},j_{m+1},\ldots,j_{n-1})\f$ giving the columns of the non-basic variables.

*/