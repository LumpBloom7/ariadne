/***************************************************************************
 *            exceptions.h
 *
 *  Copyright  2005-6  Pieter Collins, Alberto Casagrande
 *  Email  Pieter.Collins@cwi.nl, casagrande@dimi.uniud.it
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
 
/*! \file base/exceptions.h
 *  \brief Exceptions, error handling and assertions.
 */

#ifndef ARIADNE_BASE_EXCEPTIONS_H
#define ARIADNE_BASE_EXCEPTIONS_H

#include <stdexcept>
#include <iosfwd>

#include "../throw.h"
#include "base/types.h"

namespace Ariadne {
    
  //@{ \name Exceptions

  /*! \brief The method or function has not been implemented yet. */
  class NotImplemented : public std::logic_error {
   public:
    NotImplemented(const std::string& str) : std::logic_error(str) { }
  };

  /*! \brief The implementation of a virtual method has been deferred to a derived class. 
   *  However, the method is not critical for the functioning of the class, so it is not
   *  declared pure virtual. 
   */
  class DeferredImplementation : public std::logic_error {
   public:
    DeferredImplementation(const std::string& str) : std::logic_error(str) { }
  };

  /*! \brief The method or function has been deprecated. */
  class Deprecated : public std::logic_error {
   public:
    Deprecated(const std::string& str) : std::logic_error(str) { }
  };

  //@}

  namespace Base {
   
    //@{ \name Exceptions

    /*! \brief A stream or string input was invalid. */
    class InvalidInput : public std::runtime_error {
     public:
      InvalidInput(const std::string& str) : std::runtime_error(str) { }
    };
      
    /*! \brief The index to an arrayed object was invalid. */
    struct InvalidIndex : public std::runtime_error {
      InvalidIndex(const std::string& str) : std::runtime_error(str) { }
    };
      
    /*! \brief The length of an array was invalid. */
    class LengthError : public std::runtime_error {
     public:
      LengthError(const std::string& str) : std::runtime_error(str) { }
    };
      
    //@}
   
  }
  
}

#define ARIADNE_CHECK_EQUAL_ARRAY_SIZES(ary1,ary2,func) \
  { if((ary1).size()!=ary2.size()) { ARIADNE_THROW(LengthError,func,#ary1"="<<ary1<<", "#ary2"="<<ary2); } }
        
#define ARIADNE_CHECK_ARRAY_INDEX(ary,i,func) \
  { if((ary).size()<=i) { ARIADNE_THROW(InvalidIndex,func,#ary".size()="<<(ary).size()<<", "#i"="<<i); } }

#endif /* ARIADNE_BASE_EXCEPTIONS_H */