/***************************************************************************
 *            point_list.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  Email casagrande@dimi.uniud.it   Pieter.Collins@cwi.nl
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

/*! \file point_list.h
 *  \brief A list of points in Euclidean space.
 */

#ifndef _ARIADNE_POINT_LIST_H
#define _ARIADNE_POINT_LIST_H

#include <iostream>
#include <stdexcept>

#include "../declarations.h"
#include "../linear_algebra/vector.h"
#include "../linear_algebra/matrix.h"
#include "../geometry/point.h"

namespace Ariadne {
  namespace Geometry {

    template<class R>
    class PointListIterator
    {
     public:
      PointListIterator(const LinearAlgebra::Matrix<R>& A, size_type j) : _A(A), _j(j) { }
      Point<R> operator*() const { return Point<R>(column(_A,_j)); }
      void operator++() { ++_j; }
      bool operator==(const PointListIterator& other) { 
        return (this->_j==other._j) && (&(this->_A)==&(other._A)); }
      bool operator!=(const PointListIterator& other) { 
        return !(*this==other); }
     private:
      const LinearAlgebra::Matrix<R>& _A;
      size_type _j;
    };
    
    template<class R>
    class PointListReference
    {
     public:
      PointListReference(PointList<R>& pl, size_type j) : _pl(pl), _j(j) { }
      void operator=(const Point<R>& pt) {
        for(dimension_type i=0; i!=_pl.dimension(); ++i) { _pl._pts(i,_j)=pt[i]; } }
      R& operator[] (size_type i) { return _pl._pts(i,_j); }
      operator Point<R>() const { return _pl.get(_j); }
     private:
      PointList<R>& _pl;
      size_type _j;
    };
    
    /*!\brief A list of points in Euclidean space. 
     *
     * All points must have the same dimension.
     * Space is allocated for the entire list, which is more efficient than
     * allocating space for each point independently.
     *
     * A %PointList is the preferred class to hold lists of points of the same
     * dimension, such as the vertices of a polyhedron.
     * Duplicates are allowed, and in general the points are not ordered.
     *
     * A lexicographic sort() function will be made available.
     *
     * A FuzzyPointList class also exists.
     *
     * The data is stored as a (d+1)*m matrix, where d is the dimension of the
     * space and m is the number of points.
     * This is to facilitate linear programming and double description
     * method transformations. 
     * Hence the \a i th element of the \a j th point
     * is stored in the \a (i,j) th element of the matrix.
     */
    template<class R>
    class PointList
    {
      friend class PointListReference<R>;
      friend class PointListIterator<R>;
     public:
      typedef PointListIterator<R> const_iterator;
      typedef PointListIterator<R> iterator;
      
      /*!\brief Construct an empty point list, able to hold points of dimension \a d. */
      explicit PointList(dimension_type d=0) : _size(0), _pts(d+1,1) { }
      /*!\brief Construct a list consisting of \a n copies of the origin in dimension \a d. */
      explicit PointList(dimension_type d, size_type n) : _size(n), _pts(d+1,n) { }
      /*!\brief Construct from a matrix \a G whose columns contain the position vectors of the points padded by a one. */
      explicit PointList(const LinearAlgebra::Matrix<R>& G) : _size(G.number_of_columns()), _pts(G) { }
      
      /*!\brief Convert from a point list with another data type. */
      template<class Rl> PointList(const PointList<Rl>& ptl) : _size(ptl.size()), _pts(ptl.generators()) { }
      
      /*!\brief Copy constructor. */
      PointList(const PointList& ptl) : _size(ptl._size), _pts(ptl.generators()) { }
      /*!\brief Assignment operator. */
      PointList& operator=(const PointList& ptl) {
        if(this!=&ptl) { this->_size=ptl.size(); this->_pts=ptl.generators(); } 
        return *this; }
      
      /*!\brief The dimension of the points in the list. */
      dimension_type dimension() const { return _pts.number_of_rows()-1; }
      /*!\brief The number of points in the list. */
      size_type size() const { return _size; }
      /*!\brief The number of points which the list can hold without a reallocation of memory. */
      size_type capacity() const { return _pts.number_of_columns(); }
      /*!\brief Reserve space for at least \a n points. */
      void reserve(size_type n);
      
      /*!\brief A matrix expression whose columns are the position vectors of the points. */
      const LinearAlgebra::Matrix<R>& generators() const;

#ifdef DOXYGEN
      /*!\brief The \a j th point in the list. */
      Point<R> operator[] (const size_type j) const;
      /*!\brief A reference to the \a j th point in the list. */
      Point<R>& operator[] (const size_type j);
#else
      Point<R> operator[] (const size_type j) const;

      PointListReference<R> operator[] (const size_type j) { 
        return PointListReference<R>(*this,j);}
#endif

      /*!\brief Adjoin \a pt to the end of the list. */
      void push_back(const Point<R>& pt);
      /*!\brief Remove the last point from the list. */
      void pop_back() { --this->_size; }
      
      /*!\brief Adjoin \a pt to the list. */
      void adjoin(const Point<R>& pt) { this->push_back(pt); }
      /*!\brief Remove the last point from the list. */
      void clear() { this->_size=0; }

      
      /*!\brief A constant iterator to the beginning of the list. */
      const_iterator begin() const { return const_iterator(this->_pts,0); }
      /*!\brief A constant iterator to the end of the list. */
      const_iterator end() const { return const_iterator(this->_pts,this->_size); }
      
      /*!\brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
     private:
      void _set(size_type j, dimension_type i,const R& x) { _pts(i,j)=x; }
      Point<R> _get(size_type j) const { return Point<R>(column(_pts,j)); }
      R _get(size_type j, dimension_type i) const { return _pts(i,j); }

     private:
      size_type _size;
      LinearAlgebra::Matrix<R> _pts;
    };
    
    template<class R> inline
    Point<R> 
    PointList<R>::operator[] (const size_type j) const
    {
      Point<R> result(this->dimension());
      for(size_type i=0; i!=this->dimension(); ++i) {
        result[i]=this->_pts(i,j);
      }
      return result;
    }
    
    template<class R> 
    inline
    std::ostream& 
    operator<<(std::ostream& os, const PointList<R>& pts) {
      return pts.write(os);
    }
    
    
    
   
      
  }
}

#endif /* _ARIADNE_POINT_LIST_H */
