/***************************************************************************
 *            virtual_machine.code.h
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it  Pieter.Collins@cwi.nl
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
 
#include <iostream>
#include <sstream>
#include <cassert>
#include <cstring>

#include "../base/array.h"
#include "../base/stlio.h"
#include "../base/exceptions.h"
#include "../numeric/rational.h"
#include "../system/virtual_machine.h"
#include "../output/logging.h"


namespace {

using namespace Ariadne;
using namespace Ariadne::Numeric;
using Ariadne::Base::array;
using Ariadne::System::VirtualMachine;


template<class X> inline
void evaluate(const array<VirtualMachine::ByteCode>& ops, X** args, float_tag)
{
  std::vector<X> stack;
  VirtualMachine::Index index;
  int value;
  X tmp;
  for(array<VirtualMachine::ByteCode>::const_iterator op_iter=ops.begin();
      op_iter!=ops.end(); ++op_iter)
  {
    ARIADNE_LOG(7,op_iter->op);
    switch(op_iter->op) {
    case VirtualMachine::PUSH:
      index=(++op_iter)->ind;
      stack.push_back(args[index[0]][index[1]]); 
      break;
    case VirtualMachine::PULL:
      index=(++op_iter)->ind;
      args[index[0]][index[1]]=stack.back(); 
      stack.pop_back();
      break;
    case VirtualMachine::CONST:
      value=(++op_iter)->val;
      // We use the "copy and assign" idiom to prevent raw 
      // construction of differentials.
      stack.push_back(args[0][0]); 
      stack.back()=value;
      break;
    case VirtualMachine::POS:
      stack[stack.size()-1]=(+stack[stack.size()-1]);
      break;
    case VirtualMachine::NEG:
      stack[stack.size()-1]=(-stack[stack.size()-1]);
      break;
    case VirtualMachine::ADD:
      stack[stack.size()-2]=stack[stack.size()-2]+stack[stack.size()-1];
      stack.pop_back();
      break;
    case VirtualMachine::SUB:
      stack[stack.size()-2]=stack[stack.size()-2]-stack[stack.size()-1];
      stack.pop_back();
      break;
    case VirtualMachine::MUL:
      stack[stack.size()-2]=stack[stack.size()-2]*stack[stack.size()-1];
      stack.pop_back();
      break;
    case VirtualMachine::DIV:
      stack[stack.size()-2]=stack[stack.size()-2]/stack[stack.size()-1];
      stack.pop_back();
      break;
    case VirtualMachine::POW:
      //std::cerr << "POW"<<std::flush;
      value=(++op_iter)->val;
      stack[stack.size()-1]=Numeric::pow(stack[stack.size()-1],value);
      break;
    case VirtualMachine::MIN:
      stack[stack.size()-2]=Numeric::min(stack[stack.size()-2],stack[stack.size()-1]);
      stack.pop_back();
      break;
    case VirtualMachine::MAX:
      stack[stack.size()-2]=Numeric::max(stack[stack.size()-2],stack[stack.size()-1]);
      stack.pop_back();
      break;
    case VirtualMachine::ABS:
      stack[stack.size()-1]=Numeric::abs(stack[stack.size()-1]);
      break;
    case VirtualMachine::EXP:
      stack[stack.size()-1]=Numeric::exp(stack[stack.size()-1]);
      break;
    case VirtualMachine::LOG:
      stack[stack.size()-1]=Numeric::log(stack[stack.size()-1]);
      break;
    case VirtualMachine::SIN:
      stack[stack.size()-1]=Numeric::sin(stack[stack.size()-1]);
      break;
    case VirtualMachine::COS:
      stack[stack.size()-1]=Numeric::cos(stack[stack.size()-1]);
      break;
    case VirtualMachine::TAN:
      stack[stack.size()-1]=Numeric::tan(stack[stack.size()-1]);
      break;
    case VirtualMachine::ASIN:
      stack[stack.size()-1]=Numeric::asin(stack[stack.size()-1]);
      break;
    case VirtualMachine::ACOS:
      stack[stack.size()-1]=Numeric::acos(stack[stack.size()-1]);
      break;
    case VirtualMachine::ATAN:
      stack[stack.size()-1]=Numeric::atan(stack[stack.size()-1]);
      break;
    }
  }
}


template<class Q> inline
void evaluate(const array<VirtualMachine::ByteCode>& ops, Q** args, rational_tag)
{
  std::vector<Q> stack;
  VirtualMachine::Index index;
  int value;
  for(array<System::VirtualMachine::ByteCode>::const_iterator op_iter=ops.begin();
      op_iter!=ops.end(); ++op_iter)
  {
    switch(op_iter->op) {
    case VirtualMachine::PUSH:
      index=(++op_iter)->ind;
      stack.push_back(args[index[0]][index[1]]); 
      break;
    case VirtualMachine::PULL:
      index=(++op_iter)->ind;
      args[index[0]][index[1]]=stack.back(); 
      stack.pop_back();
      break;
    case VirtualMachine::POS:
      stack[stack.size()-1]=(+stack[stack.size()-1]);
      break;
    case VirtualMachine::NEG:
      stack[stack.size()-1]=(-stack[stack.size()-1]);
      break;
    case VirtualMachine::ADD:
      stack[stack.size()-2]=stack[stack.size()-2]+stack[stack.size()-1];
      stack.pop_back();
      break;
    case VirtualMachine::SUB:
      stack[stack.size()-2]=stack[stack.size()-2]-stack[stack.size()-1];
      stack.pop_back();
      break;
    case VirtualMachine::MUL:
      stack[stack.size()-2]=stack[stack.size()-2]*stack[stack.size()-1];
      stack.pop_back();
      break;
    case VirtualMachine::DIV:
      stack[stack.size()-2]=stack[stack.size()-2]/stack[stack.size()-1];
      stack.pop_back();
      break;
    case VirtualMachine::POW:
      value=(++op_iter)->val;
      stack[stack.size()-1]=Numeric::pow(stack[stack.size()-1],value);
      break;
    case VirtualMachine::MIN:
      stack[stack.size()-2]=Numeric::min(stack[stack.size()-2],stack[stack.size()-1]);
      stack.pop_back();
      break;
    case VirtualMachine::MAX:
      stack[stack.size()-2]=Numeric::max(stack[stack.size()-2],stack[stack.size()-1]);
      stack.pop_back();
      break;
    case VirtualMachine::ABS:
      stack[stack.size()-1]=Numeric::abs(stack[stack.size()-1]);
      break;
    default:
      throw std::runtime_error("evaluate: operation not permitted for rational numbers");
    }
  }
}


template<class X> inline
void 
evaluate(const array<VirtualMachine::ByteCode>& program, X** arguments) 
{
  typedef typename Numeric::traits<X>::number_type number_type;
  typedef typename Numeric::traits<number_type>::type tag;
  ::evaluate(program,arguments,tag());
}


}


namespace Ariadne {

template<class X>
void 
System::VirtualMachine::evaluate(const array<ByteCode>& program, X** arguments) const
{
  ::evaluate(program,arguments);
}

}






