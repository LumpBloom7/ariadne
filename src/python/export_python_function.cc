/***************************************************************************
 *            python/export_python_function.cc
 *
 *  Copyright  2007  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is diself_ns::stributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "python/float.h"
#include "python/read_array.h"

#include "numeric/rational.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/matrix.h"
#include "function/affine_derivative.h"
#include "function/taylor_derivative.h"
#include "function/function_interface.h"

#include <boost/python.hpp>
using namespace boost::python;

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Function;
using namespace Ariadne::Python;

template<class X> 
array<X>
extract_array(const boost::python::object& obj)
{
  list elements=extract<list>(obj);
  int n=len(elements);
  array<X> result(n);
  for(int i=0; i!=n; ++i) {
    X value=extract<X>(elements[i]);
    result[i]=value;
  }
  return result;
}

template<class X> void read(Vector<X>& v, const object& obj) {
  read_array(v.data(),obj); }

template<class X> void read(AffineDerivative<X>& ad, const object& obj) {
  array< AffineVariable<X> > ava=extract_array< AffineVariable<X> >(obj); 
  ad.resize(ava.size(),ava[0].argument_size());
  for(size_type i=0; i!=ava.size(); ++i) { ad[i]=ava[i]; } 
}

template<class X> 
void read(TaylorDerivative<X>& td, const object& obj) {
  array< TaylorVariable<X> > tva=extract_array< TaylorVariable<X> >(obj); 
  td=TaylorDerivative<X>(tva.size(),tva[0].argument_size(),tva[0].degree()); 
  for(size_type i=0; i!=tva.size(); ++i) { td[i]=tva[i]; } 
}

template<class R>
class PythonFunction
  : public FunctionInterface<R>
{
  typedef typename Numeric::traits<R>::arithmetic_type F;
 public:
  PythonFunction(std::string& nm, size_type rs, size_type as, const object& pyf) : _name(nm), _result_size(rs), _argument_size(as), _pyf(pyf) { }
  PythonFunction(size_type rs, size_type as, const object& pyf) : _name(), _result_size(rs), _argument_size(as), _pyf(pyf) { }
  PythonFunction(const object& pyf)
    : _name(), 
      _result_size(extract<int>(pyf.attr("result_size"))), 
      _argument_size(extract<int>(pyf.attr("argument_size"))), 
      _pyf(pyf) { }

  PythonFunction<R>* clone() const { return new PythonFunction<R>(*this); }
  virtual size_type result_size() const { return this->_result_size; }
  virtual size_type argument_size() const { return this->_argument_size; }
  virtual smoothness_type smoothness() const { return 255; }

  virtual Vector<F> evaluate (const Vector<F>& x) const { 
    Vector<F> r; 
    read(r,this->_pyf(x)); 
    return r; }
  virtual Matrix<F> jacobian (const Vector<F>& x) const { 
    TaylorDerivative<F> rj(this->_result_size,this->_argument_size,1u); 
    TaylorDerivative<F> aj=TaylorDerivative<F>::variable(x.size(),x.size(),1u,x); 
    read(rj,this->_pyf(aj)); 
    return rj.jacobian(); }
  virtual TaylorDerivative<F> derivative (const Vector<F>& x, const smoothness_type& d) const {  
    TaylorDerivative<F> rd(this->_result_size,this->_argument_size,d); 
    TaylorDerivative<F> ad=TaylorDerivative<F>::variable(x.size(),x.size(),d,x); 
    read(rd,this->_pyf(ad)); 
    return rd; }

  virtual std::ostream& write(std::ostream& os) const { 
    os << "Function( ";
    if(this->_name.size()>0) { os << "name=" << this->_name << ", "; }
    os << "result_size="<<this->_result_size;
    os << ", argument_size="<<this->_argument_size;
    return os << " )"; }
 private:
  std::string _name;
  size_type _result_size;
  size_type _argument_size;
  boost::python::object _pyf;
};


template<class R>
void export_python_function() 
{
  typedef typename Numeric::traits<R>::arithmetic_type F;

  class_<PythonFunction<R>, bases< FunctionInterface<R> > > python_function_class("AriadneFunction", init<object>());
  python_function_class.def(init<uint,uint,object>());
  python_function_class.def("argument_size", &PythonFunction<R>::argument_size);
  python_function_class.def("result_size", &PythonFunction<R>::result_size);
  python_function_class.def("smoothness", &PythonFunction<R>::smoothness);
  python_function_class.def("__call__",&PythonFunction<R>::operator());
  python_function_class.def("evaluate",&PythonFunction<R>::evaluate);
  python_function_class.def("jacobian",&PythonFunction<R>::jacobian);
  python_function_class.def("derivative",&PythonFunction<R>::derivative);
  python_function_class.def(self_ns::str(self));
}

template void export_python_function<Rational>();
template void export_python_function<FloatPy>();