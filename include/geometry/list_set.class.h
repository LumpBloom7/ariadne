/***************************************************************************
 *            list_set.class.h
 *
 *  Copyright  2005-7  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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
 

#ifndef ARIADNE_LIST_SET_CLASS_H
#define ARIADNE_LIST_SET_CLASS_H

#include <iosfwd>
#include <iostream>
#include <exception>
#include <stdexcept>

#include <vector>

#include "base/tribool.h"
#include "geometry/set_interface.h"

namespace Ariadne {
  namespace Geometry {

    class denotable_set_tag;
    template<class R> class Box;
    template<class T, class BS> class TimedSet;

    template<class BS> class ListSet;
    template<class T, class BS> class ListSet< TimedSet<T,BS> >;
  
    /*!\ingroup DenotableSet
     * \ingroup List
     * \brief A finite union of basic sets, represented as a sequence.
     *
     * A list set is the simplest type of denotable set class. It can hold
     * arbitrary lists of basic sets of the same type. Hence, the %ListSet class
     * also takes a template parameter which is the type of basic set contained
     * in the list.
     *
     * A list set is ordered by the order of insertion. Hence, as well as the
     * standard adjoin() method for denotable sets, a %ListSet also provides
     * the STL methods push_back() and pop_back().
     *
     * \internal ListSet is parameterised by the basic set type <class BS> instead
     * of the real type and basic set template <class R, template<class> class BS>
     * since this is more expressive; using the latter, we could not use basic set
     * types which do not take a single real parameter, and there are also times
     * when we do not know the template-id, only the basic set type.
     */
    template<class BS>
    class ListSet 
    {
     private:
      typedef typename BS::real_type R;

      /* List of basic sets. Note that std::vector provides a
       * reserve(size_type) method to increase the capacity.
       */
      dimension_type _dimension;
      std::vector< BS > _vector;

     public:
      /*! \brief A tag describing the type of set. */
      typedef denotable_set_tag set_category;
      /*!\brief The type of denotable real number used to represent points in the space. */
      typedef R real_type;
      /*!\brief The type of point contained by the set. */
      typedef Point<R> state_type;
      /*!\brief The type of basic set making up the denotable set. */
      typedef BS basic_set_type;
      /*!\brief The type of basic set in the list of sets. */
      typedef BS value_type;

      typedef typename std::vector<basic_set_type>::const_iterator const_iterator;
      typedef typename std::vector<basic_set_type>::const_iterator iterator;

     public:
      /*! \brief Destructor. */
      virtual ~ListSet();      

      /*! \brief An empty list set which can hold sets of an unspecified dimension. */
      ListSet();

      /*! \brief An empty list set which can only hold sets of dimension \a n. */
      ListSet(dimension_type n);

      /*! \brief A list set containing the element \a bs. */
      ListSet(const BS& bs);

      /*! \brief Convert from a list set of another type. */
      template<class BST> ListSet(const ListSet<BST>& ls);

      /*! \brief Convert from a pair of iterators. */
      template<class Iter> ListSet(Iter first, Iter last);



      /*! \brief Returns the number of basic sets forming this object. */
      size_type size() const;

      /*! \brief Accesses the i-th BasicSet. */
      const BS& operator[](size_type index) const;

      /*! \brief Make the set empty. */
      void clear();
      
      /*! \brief A constant iterator to the beginning of the list of basic sets. */
      const_iterator begin() const;

      /*! \brief A constant iterator to the end of the list of basic sets. */
      const_iterator end() const;

      /*! \brief Returns the denotable set's space dimension. */
      dimension_type dimension() const;

      /*! \brief Removes a set from the list and return it. */
      BS pop();

      /*! \brief Adjoins (makes union with) a basic set. */
      void adjoin(const BS& bs);

      /*! \brief Adjoins (makes union with) another list set. */
      void adjoin(const ListSet<BS>& ls);
      
      /*! \brief Adjoins (makes union with) a set. */
      template<class S> void adjoin(const S& s);

      //@{
      //! \name Input/output operators
      /*! \brief Output a summary of the set details. */
      std::string summary() const;
      /*! \brief Write to an output stream. */
      std::ostream& write(std::ostream& os) const;
      /*! \brief Read from an input stream. */
      std::istream& read(std::istream& is);
      //@}

#ifdef DOXYGEN
      //@{ 
      //! \name Geometric binary predicates
      /*! \brief Tests disjointness.
       */
      friend tribool disjoint(const ListSet<BS>& A,
                              const ListSet<BS>& B);

      /*! \brief Tests inclusion of \a A in \a B.
       */
      friend tribool subset(const ListSet<BS>& A,
                            const ListSet<BS>& B);
      //@}
      
      
      //@{
      //! \name Geometric binary operations
      /*! \brief The union of \a A and \a B.
       * <br><br>
       * Note that 'union' is a reserved word in C++.
       */
      //FIXME: Compiler doesn't like this
      //friend ListSet<BS> join<> (const ListSet<BS>& A,
      //                           const ListSet<BS>& B);
      friend ListSet<BS> join(const ListSet<BS>& A,
                              const ListSet<BS>& B);

      /*! \brief The closure of intersection of the interior of \a A with the interior of \a B.
       */
      friend ListSet<BS> open_intersection(const ListSet<BS>& A,
                                           const ListSet<BS>& B);
      /*! \brief The union of all cells of \a A which can be proved to be subsets of \a B.
       */
      friend ListSet<BS> inner_intersection(const ListSet<BS>& A,
                                            const SetInterface<R>& B);
      /*! \brief The union of all cells of \a A which can be proved to intersect \a B.
       */
      friend ListSet<BS> lower_intersection(const ListSet<BS>& A,
                                            const SetInterface<R>& B);

      /*! \brief The union of all cells of \a A which have not been proved to be disjoint from \a B.
       */
      friend ListSet<BS> outer_intersection(const ListSet<BS>& A,
                                            const SetInterface<R>& B);
      //@}
#endif      
     private:
      template<class S> void adjoin(const S& bs, basic_set_tag);
      template<class DS> void adjoin(const DS& ds, denotable_set_tag);

     private:
      static void _instantiate();
    };

  
    template<class BS>
    std::ostream& operator<<(std::ostream& os, const ListSet<BS>& A);


    template<class BS>
    std::istream& operator>>(std::istream& is, ListSet<BS>& A);

    template<class BS>
    Box<typename BS::real_type> bounding_box(const ListSet<BS>& A);

    template<class BS, class R>
    tribool disjoint(const ListSet<BS>& A, const Box<R>& B);

    template<class BS>
    ListSet<BS> open_intersection(const ListSet<BS>& A, const ListSet<BS>& B);
  
    template<class BS1, class BS2>
    tribool disjoint(const ListSet<BS1>& A, const ListSet<BS2>& B);

    template<class BS1, class BS2>
    tribool subset(const ListSet< BS1 >& A, const BS2& B);
    
    template<class BS>
    ListSet<BS> join(const ListSet<BS>& A, const ListSet<BS>& B);
    

    template<class BS>
    ListSet<BS> open_intersection(const ListSet<BS>& A, const ListSet<BS>& B);
  
  
  }
}

#endif /* ARIADNE_LIST_SET_CLASS_H */