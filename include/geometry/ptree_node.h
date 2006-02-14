/***************************************************************************
 *            ptree_node.h
 *
 *  18 January 2005
 *  Copyright  2004,2005  Alberto Casagrande, Pieter Collins
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


#ifndef _ARIADNE_PTREE_NODE_H
#define _ARIADNE_PTREE_NODE_H
#include <bitset>
#include <vector>
#include <stdexcept>

#include <iostream>
#include <fstream>  
#include <iomanip> 
#include <string>

#include <geometry/rectangle.h>

#define MAX_DIM 12
#define MAX_BITSET 1<<MAX_DIM

namespace Ariadne {
namespace Geometry {

template <typename R> class PTree;
class PTreeNode;

std::ostream& operator<<(std::ostream &os, const PTreeNode &node);


/*! \brief Intersect two subtrees. */
PTreeNode intersect(const PTreeNode &A, const PTreeNode &B);

/*! \brief Intersect two subtrees. */
PTreeNode join(const PTreeNode &A, const PTreeNode &B);

/*! \brief A node of a PTree partition tree. */
class PTreeNode {
  typedef std::bitset<MAX_BITSET> space_type;
 protected:
  /*! \brief  The flag for leafs identification. */
  bool _leaf;
  
  /*! \brief  This flag is true if the node represents an empty 
   * space. */
  bool _empty;
  
  /*! \brief  The node's depth. */
  size_t _depth;
  
  /*! \brief  The node's level. */
  size_t _level;
  
  /*! \brief  The space dimension. */
  size_t _dimension;
  
  /*! \brief  The 0/1 representation of the space. */
  space_type _space;
  
  /*! \brief  The list of not empty node's sons. 
   *
   * Notice that |_sons| is equal to _space.count() if
   * _leaf=false. Otherwise, if _leaf=true, it is equal to 0.
   */
  std::vector<PTreeNode *> _sons;
  
  /*! \brief  Returns true if the subtree is empty. */
  inline bool empty() const {
    return (this->_empty);
  }
  
  inline void set_not_empty() {
    this->_empty=false;
    this->_leaf=false;
  }
  
  /*! \brief  Returns true if the subtree is full. */
  inline bool full() const {
    return ((!this->_empty)&&(this->_leaf));
  }
  
  inline void set_empty() {
    this->_empty=true;
    (this->_space).reset();
    this->_leaf=true;
  }
  
  inline void set_full() {
    
    this->_clear_sons();
    this->_level=0;
    
    this->_empty=false;
    (this->_space).reset();
    this->_leaf=true;
  }
  
  /*! \brief  Returns the number of bits used by _space. */
  inline size_t _bits_per_space_member() const {
    return (1<<(this->_dimension));
  }
  
  /*! \brief Increases node's depth. */
  inline void increase_depth() {
    
    this->_depth++;
    
    for (size_t i=0;i<(this->_sons).size();i++) {
      if (this->_sons[i]!=NULL) {
        (this->_sons[i])->increase_depth();
      }
    }
    
  }
  
  /*! \brief  Evaluates node's level using son levels. */
  inline size_t eval_level_by_sons() {
    
#ifdef DEBUG
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
#endif
    
    this->_level=0;
    
    /* compare this->_level with each son level + 1*/
    for (size_t i=0;i<(this->_sons).size();i++) {
      if (this->_sons[i]!=NULL) {
        if (((this->_sons[i])->_level+1)>this->_level) {
          
          /* new this->_level is son level + 1*/
          this->_level=(this->_sons[i])->_level+1;
        }
      }
    }
    
#ifdef DEBUG
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
#endif
    
    return this->_level;
  }
  
  inline bool _check_emptyness_and_set_leaf(){
    
    if (this->_space.any()) {
      this->set_not_empty();
    }
    
    return this->_leaf;
  }
  
  /*! \brief  Intersects node's sons. */
  inline void intersect_sons(const PTreeNode &A, const PTreeNode &B) {
    
    
    /* intersect used quadrant */
    this->_space=A._space & B._space;
    
    if (this->_check_emptyness_and_set_leaf()) {
      return;
    }
    
    /* for each quadrant i*/
    for (size_t i=0;i<A._bits_per_space_member();i++) {
      
      /* if i is used by A and by B */
      if ((A._space).test(i)&&(B._space).test(i)) {
        
        /* intesect the two quardants */
        PTreeNode *C=new PTreeNode(intersect(*A._sons[i], *B._sons[i]));
        
        /* if the intersection is not empty */
        if (!C->empty()) {
          
          /* insert the intersection at the end of son 
           * vector */
          this->_sons[i]=C;
        } else {
          
          /* if the intersection 
           * is empty, reset the 
           * i quadrant into 
           * _space */
          this->_sons[i]=NULL;
          (this->_space).reset(i);
          
          delete C;
        }
        
      }
      
    }
    
    /* if this->_space is empty this is an empty leaf */
    if ((this->_space).none()) {
      this->set_empty();
      
      return;
    }
    
  }
  
  inline void compact_node() {
    
    /* if every son contains something */
    if ((this->_space).count() == this->_bits_per_space_member()) {
      
      /* check if every son is a full leaf */
      bool full_leaf=true;
      
      for (size_t i=0;i<this->_bits_per_space_member();i++) {
        full_leaf= full_leaf && (this->_sons[i])->full();
      }
      
      /* if every son is a full leaf */
      if (full_leaf==true) {
        
        /* transform this node into a full leaf*/
        this->set_full();
      }
    }
  }
  
  inline void inplace_join_sons(const PTreeNode &A) {
    
    /* for each quadrant i*/
    for (size_t i=0;i<A._bits_per_space_member();i++) {
      
      /* if i is used by A */
      if ((A._space).test(i)) {
        
        /* if i is used by this object */
        if (this->_space.test(i)) {
          
          /* join the two quadrants */
          (this->_sons[i])->
            inplace_union(*A._sons[i]);
          
        } else {
          /* if i is not used by this 
           * object */
          
          /* add a copy of the i quadrant of A at the end of
           * the son vector */
          this->_sons[i]=new PTreeNode(*A._sons[i]);
        }
      }
      
    }
    
    /* join used quadrant */
    this->_space=A._space | this->_space;
    
    if (this->_check_emptyness_and_set_leaf()) {
      return;
    }
    
    /* compact node */
    this->compact_node();
    
  }
  
  
  inline void join_sons(const PTreeNode &A, const PTreeNode &B) {
    
    /* join used quadrant */
    this->_space=A._space | B._space;
    
    if (this->_check_emptyness_and_set_leaf()) {
      return;
    }
    
    /* for each quadrant i*/
    for (size_t i=0;i<A._bits_per_space_member();i++) {
      
      /* if i is used by A */
      if ((A._space).test(i)) {
        
        /* if i is used by B */
        if ((B._space).test(i)) {
          
          /* join the two quadrants */
          
          this->_sons[i]=new PTreeNode(join(*A._sons[i], *B._sons[i]));
          
        } else {
          /* if i is not used by B */
          
          /* add a copy of the i quadrant of A at the end of
           * the son vector */
          this->_sons[i]=new PTreeNode(*A._sons[i]);
        }
        
      } else {
        /* if i is not used by A */
        
        /* if i is used by B */
        if ((B._space).test(i)) {
          
          /* add a copy of the i quadrant of B at the end of
           * the son vector */
          this->_sons[i]=new PTreeNode(*B._sons[i]);
          
        } else {
          
          this->_sons[i]=NULL;
        }
      }
      
    }
    
    /* compact node */
    this->compact_node();
    
  }
  
  void _clear_sons() {
    
    /* delete sons */
    for (size_t i=0;i<(this->_sons).size();i++) {
      delete (this->_sons[i]);
    }
    
    (this->_sons).clear();
  }
  
 public:
  /*! \brief A class costructor. */
  PTreeNode(): _leaf(true), _empty(true), _depth(0), 
                       _level(0), _dimension(0) {
    
#ifdef DEBUG
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
#endif
    
    (this->_space).reset();
    (this->_sons).clear();
    
    if (this->_dimension>MAX_DIM) 
      throw std::out_of_range("Wrong space dimension"); 
    
#ifdef DEBUG
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
#endif
  }
  
  /*! \brief A class costructor. */
  PTreeNode(const size_t dimension): _leaf(true), _empty(true), 
                                             _depth(0), _level(0),_dimension(dimension) {
    
#ifdef DEBUG
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
#endif
    
    if (this->_dimension>MAX_DIM) 
      throw std::out_of_range("Wrong space dimension"); 
    
    
    (this->_space).reset();
    (this->_sons).resize(this->_bits_per_space_member());
    
    for (size_t i=0;i<(this->_sons).size();i++) {
      this->_sons[i]=NULL;
    }
    
#ifdef DEBUG
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
#endif
    
  }
  
  /*! \brief A class costructor. */
  PTreeNode(const PTreeNode &A): _leaf(A._leaf), 
                                                 _empty(A._empty), _depth(A._depth),  
                                                 _level(A._level), 
                                                 _dimension(A._dimension), _space(A._space) {
    
#ifdef DEBUG
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
#endif
    
    (this->_sons).resize(A._bits_per_space_member());
    
    for (size_t i=0;i<(A._sons).size();i++) {
      
      if (A._sons[i]!=NULL) {
        this->_sons[i]=new PTreeNode(*A._sons[i]);
      } else {
        this->_sons[i]=NULL;
      }
    }
    
#ifdef DEBUG
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
#endif
    
  }
  
  /*! \brief A class destructor. */
  ~PTreeNode(){
    
    this->_clear_sons();
  }
  
  inline PTreeNode &operator=(const PTreeNode &A){
    
#ifdef DEBUG
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
#endif
    
    this->_clear_sons();
    
    this->_leaf=A._leaf;
    this->_depth=A._depth;
    this->_level=A._level;
    this->_space=A._space;
    this->_dimension=A._dimension;
    this->_empty=A._empty;
    
    (this->_sons).resize(A._bits_per_space_member());
    
    for (size_t i=0;i<(A._sons).size();i++) {
      if (A._sons[i]!=NULL) {
        this->_sons[i]=new PTreeNode(*A._sons[i]);
      } else {
        this->_sons[i]=NULL;
      }
    }
    
#ifdef DEBUG
    std::cout << __FILE__ << ":" << __LINE__ << std::endl;
#endif
    
    return *this;
  }
  
  inline size_t dimension() const{
    return this->_dimension;
  }
  
  inline size_t dim() const{
    return this->_dimension;
  }
  
  inline PTreeNode add_upper_level(unsigned int i) {
    
    if (this->_bits_per_space_member()<=i)
      throw std::out_of_range("Wrong index"); 
    
    PTreeNode A(this->_dimension);
    
    (A._space).set(i);
    A._sons[i]=new PTreeNode(*this);
    
    A._level= this->_level+1;
    A._sons[i]->increase_depth();
    
    A.set_not_empty();
    
    return A;
  }
  
  inline PTreeNode add_lower_level(PTreeNode &A, unsigned int i) const{
    
    if (this->_dimension!=A._dimension) 
      throw std::domain_error("The object and parameter have different space dimensions");
    
    if (A._bits_per_space_member()<=i)
      throw std::out_of_range("Wrong index");
    
    PTreeNode new_node = A.add_upper_level(i);
    
    new_node= join(*this,new_node);
    
    return new_node;
  }
  
  /*! \brief Check inclusion. */
  inline bool does_include(const PTreeNode &A) const{
    
    /* intersect used quadrant */
    space_type c_space=this->_space & A._space;
    
    /* if every quadrant used by A is also used by this 
     * object */
    if ((A._space).count()==c_space.count()) {
      
      /* for each quadrant i*/
      for (size_t i=0;i<A._bits_per_space_member();i++) {
        
        /* if i is used by A */
        if ((A._space).test(i)) {
          
          if (!((this->_sons[i])->
                does_include(*A._sons[i]))) {
            
            return false;
            
          }
          
        }
        
      }
      
      /* if every subtree is included */
      return true;
    }
    
    /* if A uses a quadrant which this object does not use*/
    return false;
  }
  
  /*! \brief Inplace join. */
  void inplace_union(const PTreeNode &A){
    
    if (this->_dimension!=A._dimension) 
      throw std::domain_error("The two parameters have different space dimensions.");
    
    
    if (A.empty()||this->full()) {
      return;
    }
    
    if (this->empty()||A.full()) {
      *this=A;
      return;
    }
    
    this->inplace_join_sons(A);
    this->eval_level_by_sons();
  }
  
  template <typename R> friend class PTree;
  
  friend std::ostream& Ariadne::Geometry::operator<<(std::ostream &os, 
                                                     const PTreeNode &node);
  
  /*! \brief Intersect two subtrees. */
  friend PTreeNode intersect(const PTreeNode &A, const PTreeNode &B);
  
  /*! \brief Intersect two subtrees. */
  friend PTreeNode join(const PTreeNode &A, const PTreeNode &B);
  
  //friend class Ariadne::Geometry::IO_Operators::Exporter;
};
  
/*! \brief Intersect two subtrees. */
PTreeNode intersect(const PTreeNode &A, const PTreeNode &B) {
  
  if (A._dimension!=B._dimension) 
    throw std::domain_error("This object and parameter have different space dimensions.");
  
  
  if (A._leaf||B.empty()) {
    return B;
  } else {
    if (B._leaf||A.empty()) {
      return A;
    } else {
      
      PTreeNode C(A._dimension);
      
      C.intersect_sons(A,B);
      
      C._depth=A._depth;
      C.eval_level_by_sons();
      
      return C;
      
    }
    
  }
  
}
  
/*! \brief Join two subtrees. */
PTreeNode join(const PTreeNode &A, const PTreeNode &B){
  if (A._dimension!=B._dimension) 
    throw std::domain_error("The two parameters have different space dimensions.");
  
  if (A.empty()||B.full()) {
    return B;
  }
  else {
    if (B.empty()||A.full()) {
      return A;
    }
    else {
      PTreeNode C(A._dimension);
      
      C.join_sons(A,B);
      C._depth=A._depth;
      C.eval_level_by_sons();
      
      return C;
    }
  }
}
  
std::ostream& operator<<(std::ostream &os, const Ariadne::Geometry::PTreeNode &node) {
  os << "[ (level=" << node._level << ")(depth=" << node._depth << ")";
  
  if (node.full()) {
    os << "1]";
    return os;
  }
  
  if (node.empty()) {
    os << "0]";
    
    return os;
  }
  
  size_t idx_this=0;
    
  /* for each quadrant i*/
  for (size_t i=0;i<node._bits_per_space_member();i++) {
    /* if i is used by A */
    if ((node._space).test(i)) {
      os << *(node._sons[idx_this]);
      idx_this++;
    }
    else {
      /* if i is not used by A */
      os << "[]";
    }
  }
  os << " ]"; 
  
  return os;
}
 
}
}

#endif //_ARIADNE_PTREE_NODE_H
