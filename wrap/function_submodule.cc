/***************************************************************************
 *            function_submodule.cc
 *
 *  Copyright 2008  Pieter Collins
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

#include "utilities.h"

#include <boost/python.hpp>

#include <iostream>
#include <iomanip>

#include "array.h"
#include "container.h"
#include "numeric.h"
#include "vector.h"
#include "expansion.h"
#include "multi_index.h"
#include "taylor_model.h"
#include "differential.h"
#include "polynomial.h"
#include "affine.h"
#include "taylor_function.h"
#include "constraint.h"
#include "function.h"
#include "function_mixin.h"
#include "expression.h"
#include "space.h"
#include "assignment.h"
#include "formula.h"
#include "constraint_solver.h"

#include "../src/function_mixin.tcc"

using namespace boost::python;

namespace Ariadne {

template<>
struct from_python< MultiIndex >
{
    from_python() { converter::registry::push_back(&convertible,&construct,type_id<MultiIndex>()); }
    static void* convertible(PyObject* obj_ptr) { if (!PyTuple_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        boost::python::tuple tup=extract<boost::python::tuple>(obj_ptr);
        void* storage = ((converter::rvalue_from_python_storage<MultiIndex>*)   data)->storage.bytes;
        MultiIndex res(len(tup));
        for(uint i=0; i!=res.size(); ++i) { res.set(i,extract<uint>(tup[i])); }
        new (storage) MultiIndex(res);
        data->convertible = storage;
    }
};




template<>
struct from_python<EffectiveVectorFunction>
{
    from_python() { converter::registry::push_back(&convertible,&construct,type_id<EffectiveVectorFunction>()); }
    static void* convertible(PyObject* obj_ptr) { if (!PyList_Check(obj_ptr)) { return 0; } return obj_ptr; }
    static void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data) {
        list lst=extract<list>(obj_ptr);
        void* storage = ((converter::rvalue_from_python_storage< EffectiveVectorFunction >*)   data)->storage.bytes;
        EffectiveVectorFunction res(len(lst),0);
        for(uint i=0; i!=res.result_size(); ++i) { res.set(i,extract<EffectiveScalarFunction>(lst[i])); }
        new (storage) EffectiveVectorFunction(res);
        data->convertible = storage;
    }
};




template<class X, class D>
inline Matrix<X> get_jacobian(const Vector<D>& d) {
    const uint rs=d.size(); const uint as=d[0].argument_size();
    Matrix<X> J(rs,as);
    for(uint i=0; i!=rs; ++i) {
        for(uint j=0; j!=as; ++j) {
            J[i][j]=d[i][j];
        }
    }
    return J;
}

template<class X> std::ostream& operator<<(std::ostream& os, const Representation< ScalarFunction<X> >& frepr) {
    static_cast<const ScalarFunctionInterface<X>&>(frepr.reference()).repr(os); return os;
}

template<class X> std::ostream& operator<<(std::ostream& os, const Representation< VectorFunction<X> >& frepr) {
    static_cast<const VectorFunctionInterface<X>&>(frepr.reference()).repr(os); return os;
}

class ScalarPythonFunction
    : public ScalarFunctionMixin<ScalarPythonFunction,Real>
{
    friend class ScalarFunctionMixin<ScalarPythonFunction,Real>;
    template<class T> void _compute(T& r, const Vector<T>& a) const {
        r=boost::python::extract<T>(this->_pyf(a)); }
  public:
    ScalarPythonFunction(std::string& nm, uint as, const object& pyf) : _name(nm), _argument_size(as), _pyf(pyf) { }
    ScalarPythonFunction(uint as, const object& pyf) : _name(),  _argument_size(as), _pyf(pyf) { }
    ScalarPythonFunction(const object& pyf)
        : _name(),
          _argument_size(extract<int>(pyf.attr("argument_size"))),
          _pyf(pyf) { }

    ScalarPythonFunction* clone() const { return new ScalarPythonFunction(*this); }
    virtual uint argument_size() const { return this->_argument_size; }

    virtual Vector<Float> gradient(const Vector<Float>& x) const {
        return this->evaluate(Differential<Float>::variables(1u,x)).gradient(); }
    virtual Vector<Interval> gradient(const Vector<Interval>& x) const {
        return this->evaluate(Differential<Interval>::variables(1u,x)).gradient(); }

    virtual EffectiveScalarFunctionInterface* _derivative (uint j) const {
        ARIADNE_FAIL_MSG("Cannot get a component of a Python function"); }
    virtual std::ostream& repr(std::ostream& os) const { return os; }
    virtual std::ostream& write(std::ostream& os) const {
        os << "ScalarUserFunction( ";
        if(this->_name.size()>0) { os << "name=" << this->_name << ", "; }
        os << "argument_size="<<this->_argument_size;
        return os << " )"; }
    EffectiveScalarFunction derivative(uint j) const { return this->_derivative(j); }
  private:
    std::string _name;
    uint _argument_size;
    boost::python::object _pyf;
};


class VectorPythonFunction
    : public VectorFunctionMixin<VectorPythonFunction,Real>
{
    friend class VectorFunctionMixin<VectorPythonFunction,Real>;
    template<class T> void _compute(Vector<T>& r, const Vector<T>& a) const {
        r=boost::python::extract< Vector<T> >(this->_pyf(a)); }
  public:
    VectorPythonFunction(std::string& nm, uint rs, uint as, const object& pyf) : _name(nm), _result_size(rs), _argument_size(as), _pyf(pyf) { }
    VectorPythonFunction(uint rs, uint as, const object& pyf) : _name(), _result_size(rs), _argument_size(as), _pyf(pyf) { }
    VectorPythonFunction(const object& pyf)
        : _name(),
          _result_size(extract<int>(pyf.attr("result_size"))),
          _argument_size(extract<int>(pyf.attr("argument_size"))),
          _pyf(pyf) { }

    VectorPythonFunction* clone() const { return new VectorPythonFunction(*this); }
    virtual uint result_size() const { return this->_result_size; }
    virtual uint argument_size() const { return this->_argument_size; }

    virtual Matrix<Float> jacobian (const Vector<Float>& x) const {
        return this->evaluate(Differential<Float>::variables(1u,x)).jacobian(); }
    virtual Matrix<Interval> jacobian (const Vector<Interval>& x) const {
        return this->evaluate(Differential<Interval>::variables(1u,x)).jacobian(); }

    virtual EffectiveScalarFunctionInterface* _get(uint i) const {
        ARIADNE_FAIL_MSG("Cannot get a component of a Python function"); }
    virtual EffectiveScalarFunction operator[](uint i) const {
        ARIADNE_FAIL_MSG("Cannot get a component of a Python function"); }

    virtual std::ostream& write(std::ostream& os) const {
        os << "VectorUserFunction( ";
        if(this->_name.size()>0) { os << "name=" << this->_name << ", "; }
        os << "result_size="<<this->_result_size;
        os << ", argument_size="<<this->_argument_size;
        return os << " )"; }
  private:
    std::string _name;
    uint _result_size;
    uint _argument_size;
    boost::python::object _pyf;
};



}

using namespace Ariadne;

typedef Float F;
typedef Interval I;
typedef Vector<Float> FV;
typedef Vector<Interval> IV;
typedef Matrix<Float> FMx;
typedef Matrix<Interval> IMx;
typedef Vector< Differential<Float> > FSDV;
typedef Vector< Differential<Interval> > ISDV;
typedef Vector<IntervalTaylorModel> TMV;
typedef VectorTaylorFunction TFM;
typedef IntervalTaylorModel TM;





void export_multi_index()
{
    class_< MultiIndex > multi_index_class("MultiIndex", init<uint>());
    multi_index_class.def(init<MultiIndex>());
    multi_index_class.def("__getitem__",&MultiIndex::get);
    multi_index_class.def("__setitem__",&MultiIndex::set);
    multi_index_class.def("degree",&MultiIndex::degree);
    multi_index_class.def(self_ns::str(self));
    multi_index_class.def(self_ns::repr(self));

    from_python<MultiIndex>();
    to_python< List<MultiIndex> >();
}

template<class X>
void export_monomial()
{
    typedef ExpansionValue<X> M;
    class_< M > monomial_class(python_name<X>("Monomial"), init<MultiIndex,X>());
    monomial_class.def("key",(const MultiIndex&(M::*)()const)&M::key,return_value_policy<copy_const_reference>());
    monomial_class.def("data",(const X&(M::*)()const) &M::data,return_value_policy<copy_const_reference>());
    monomial_class.def(self_ns::str(self));
}

template<class X>
void export_polynomial()
{
    X real;

    class_< Polynomial<X> > polynomial_class(python_name<X>("Polynomial"), init< Polynomial<X> >());
    polynomial_class.def(init<uint>());
    polynomial_class.def("constant", (Polynomial<X>(*)(uint,double)) &Polynomial<X>::constant);
    polynomial_class.staticmethod("constant");
    polynomial_class.def("variable", (Polynomial<X>(*)(uint,uint)) &Polynomial<X>::variable);
    polynomial_class.staticmethod("variable");
    polynomial_class.def("coordinate", (Polynomial<X>(*)(uint,uint)) &Polynomial<X>::variable);
    polynomial_class.staticmethod("coordinate");
    polynomial_class.def("variables", (Vector< Polynomial<X> >(*)(uint)) &Polynomial<X>::variables);
    polynomial_class.staticmethod("variables");

    polynomial_class.def("argument_size", &Polynomial<X>::argument_size);
    polynomial_class.def("insert", &Polynomial<X>::insert);
    polynomial_class.def(+self);
    polynomial_class.def(-self);
    polynomial_class.def(self+self);
    polynomial_class.def(self-self);
    polynomial_class.def(self*self);
    polynomial_class.def(self+real);
    polynomial_class.def(self-real);
    polynomial_class.def(self*real);
    polynomial_class.def(self/real);
    polynomial_class.def(real+self);
    polynomial_class.def(real-self);
    polynomial_class.def(real*self);
    polynomial_class.def("__iter__",boost::python::iterator< Polynomial<X> >());
    polynomial_class.def(self_ns::str(self));
    //polynomial_class.def(self_ns::repr(self));

    to_python< Vector< Polynomial<X> > >();
}

void export_scalar_function()
{
    class_<EffectiveScalarFunction>
        scalar_function_class("EffectiveScalarFunction", init<EffectiveScalarFunction>());
    scalar_function_class.def(init<uint>());
    scalar_function_class.def("argument_size", &EffectiveScalarFunction::argument_size);
    scalar_function_class.def("derivative", &EffectiveScalarFunction::derivative);
    scalar_function_class.def("__call__", (Interval(EffectiveScalarFunction::*)(const Vector<Interval>&)const)&EffectiveScalarFunction::operator() );
    scalar_function_class.def("__call__", (Float(EffectiveScalarFunction::*)(const Vector<Float>&)const)&EffectiveScalarFunction::operator() );
    scalar_function_class.def("__call__", (IntervalDifferential(EffectiveScalarFunction::*)(const Vector<IntervalDifferential>&)const)&EffectiveScalarFunction::evaluate );
    scalar_function_class.def("__call__", (FloatDifferential(EffectiveScalarFunction::*)(const Vector<FloatDifferential>&)const)&EffectiveScalarFunction::evaluate );
    scalar_function_class.def("gradient", (Vector<Interval>(EffectiveScalarFunction::*)(const Vector<Interval>&)const)&EffectiveScalarFunction::gradient );
    scalar_function_class.def("gradient", (Vector<Interval>(EffectiveScalarFunction::*)(const Vector<Interval>&)const)&EffectiveScalarFunction::gradient );
    scalar_function_class.def("differential", (Differential<Interval>(EffectiveScalarFunction::*)(const Vector<Interval>&,Nat)const) &EffectiveScalarFunction::differential);
    scalar_function_class.def("differential", (Differential<Float>(EffectiveScalarFunction::*)(const Vector<Float>&,Nat)const) &EffectiveScalarFunction::differential);
    scalar_function_class.def("__pos__", &__pos__<EffectiveScalarFunction,EffectiveScalarFunction>);
    scalar_function_class.def("__neg__", &__neg__<EffectiveScalarFunction,EffectiveScalarFunction>);
    scalar_function_class.def("__add__", &__add__<EffectiveScalarFunction,EffectiveScalarFunction,EffectiveScalarFunction>);
    scalar_function_class.def("__sub__", &__sub__<EffectiveScalarFunction,EffectiveScalarFunction,EffectiveScalarFunction>);
    scalar_function_class.def("__mul__", &__mul__<EffectiveScalarFunction,EffectiveScalarFunction,EffectiveScalarFunction>);
    scalar_function_class.def("__div__", &__div__<EffectiveScalarFunction,EffectiveScalarFunction,EffectiveScalarFunction>);
    scalar_function_class.def("__add__", &__add__<EffectiveScalarFunction,EffectiveScalarFunction,Real>);
    scalar_function_class.def("__sub__", &__sub__<EffectiveScalarFunction,EffectiveScalarFunction,Real>);
    scalar_function_class.def("__mul__", &__mul__<EffectiveScalarFunction,EffectiveScalarFunction,Real>);
    scalar_function_class.def("__div__", &__div__<EffectiveScalarFunction,EffectiveScalarFunction,Real>);
    scalar_function_class.def("__radd__", &__radd__<EffectiveScalarFunction,EffectiveScalarFunction,Real>);
    scalar_function_class.def("__rsub__", &__rsub__<EffectiveScalarFunction,EffectiveScalarFunction,Real>);
    scalar_function_class.def("__rmul__", &__rmul__<EffectiveScalarFunction,EffectiveScalarFunction,Real>);
    scalar_function_class.def("__rdiv__", &__rdiv__<EffectiveScalarFunction,EffectiveScalarFunction,Real>);
    scalar_function_class.def("__eq__", &__eq__<RealConstraint,EffectiveScalarFunction,Real>);
    scalar_function_class.def("__le__", &__le__<RealConstraint,EffectiveScalarFunction,Real>);
    scalar_function_class.def("__ge__", &__ge__<RealConstraint,EffectiveScalarFunction,Real>);
    scalar_function_class.def("__str__", &__cstr__<EffectiveScalarFunction>);
    scalar_function_class.def("__repr__", &__crepr__<EffectiveScalarFunction>);

    scalar_function_class.def("constant", (EffectiveScalarFunction(*)(uint,Real)) &EffectiveScalarFunction::constant);
    scalar_function_class.def("coordinate", (EffectiveScalarFunction(*)(uint,uint)) &EffectiveScalarFunction::coordinate);
    scalar_function_class.staticmethod("constant");
    scalar_function_class.staticmethod("coordinate");

    def("evaluate", (Float(*)(const EffectiveScalarFunction&,const Vector<Float>&)) &evaluate);
    def("evaluate", (Interval(*)(const EffectiveScalarFunction&,const Vector<Interval>&)) &evaluate);

    def("derivative", (EffectiveScalarFunction(EffectiveScalarFunction::*)(uint)const) &EffectiveScalarFunction::derivative);

    def("pow", (EffectiveScalarFunction(*)(const EffectiveScalarFunction&,int)) &pow);
    def("rec", (EffectiveScalarFunction(*)(const EffectiveScalarFunction&)) &rec);
    def("sqr", (EffectiveScalarFunction(*)(const EffectiveScalarFunction&)) &sqr);
    def("sqrt", (EffectiveScalarFunction(*)(const EffectiveScalarFunction&)) &sqrt);
    def("exp", (EffectiveScalarFunction(*)(const EffectiveScalarFunction&)) &exp);
    def("log", (EffectiveScalarFunction(*)(const EffectiveScalarFunction&)) &log);
    def("sin", (EffectiveScalarFunction(*)(const EffectiveScalarFunction&)) &sin);
    def("cos", (EffectiveScalarFunction(*)(const EffectiveScalarFunction&)) &cos);
    def("tan", (EffectiveScalarFunction(*)(const EffectiveScalarFunction&)) &tan);

    def("lie_derivative", (EffectiveScalarFunction(*)(const EffectiveScalarFunction&,const EffectiveVectorFunction&)) &lie_derivative);

    class_<ValidatedScalarFunction> interval_scalar_function_class("ValidatedScalarFunction", init<ValidatedScalarFunction>());
    interval_scalar_function_class.def(init<EffectiveScalarFunction>());
    interval_scalar_function_class.def(init<uint>());
    interval_scalar_function_class.def("argument_size", &ValidatedScalarFunction::argument_size);
    interval_scalar_function_class.def("__str__", &__cstr__<ValidatedScalarFunction>);
    interval_scalar_function_class.def("__repr__", &__crepr__<ValidatedScalarFunction>);

    implicitly_convertible<EffectiveScalarFunction,ValidatedScalarFunction>();
}


void export_vector_function()
{

    class_<EffectiveVectorFunction>
        vector_function_class("EffectiveVectorFunction", init<EffectiveVectorFunction>());
    vector_function_class.def(init<uint,uint>());

    vector_function_class.def("result_size", &EffectiveVectorFunction::result_size);
    vector_function_class.def("argument_size", &EffectiveVectorFunction::argument_size);
    vector_function_class.def("__getitem__", &EffectiveVectorFunction::get);
    vector_function_class.def("__setitem__", &EffectiveVectorFunction::set);
    vector_function_class.def("__call__", (Vector<Interval>(EffectiveVectorFunction::*)(const Vector<Interval>&)const)&EffectiveVectorFunction::operator() );
    vector_function_class.def("__call__", (Vector<Float>(EffectiveVectorFunction::*)(const Vector<Float>&)const)&EffectiveVectorFunction::operator() );
    vector_function_class.def("__call__", (Vector<IntervalDifferential>(EffectiveVectorFunction::*)(const Vector<IntervalDifferential>&)const)&EffectiveVectorFunction::evaluate );
    vector_function_class.def("__call__", (Vector<FloatDifferential>(EffectiveVectorFunction::*)(const Vector<FloatDifferential>&)const)&EffectiveVectorFunction::evaluate );
    vector_function_class.def("jacobian", (Matrix<Interval>(EffectiveVectorFunction::*)(const Vector<Interval>&)const) &EffectiveVectorFunction::jacobian);
    vector_function_class.def("jacobian", (Matrix<Float>(EffectiveVectorFunction::*)(const Vector<Float>&)const) &EffectiveVectorFunction::jacobian);
    vector_function_class.def("differentials", (Vector<Differential<Interval> >(EffectiveVectorFunction::*)(const Vector<Interval>&,Nat)const) &EffectiveVectorFunction::differentials);
    vector_function_class.def("differentials", (Vector<Differential<Float> >(EffectiveVectorFunction::*)(const Vector<Float>&,Nat)const) &EffectiveVectorFunction::differentials);
    vector_function_class.def("__str__", &__cstr__<EffectiveVectorFunction>);
    vector_function_class.def("__repr__", &__crepr__<EffectiveVectorFunction>);

    vector_function_class.def("identity", (EffectiveVectorFunction(*)(uint)) &EffectiveVectorFunction::identity);
    vector_function_class.staticmethod("identity");

    def("evaluate", (Vector<Float>(*)(const EffectiveVectorFunction&,const Vector<Float>&)) &evaluate);
    def("evaluate", (Vector<Interval>(*)(const EffectiveVectorFunction&,const Vector<Interval>&)) &evaluate);

    def("join", (EffectiveVectorFunction(*)(const EffectiveScalarFunction&, const EffectiveScalarFunction&)) &join);
    def("join", (EffectiveVectorFunction(*)(const EffectiveVectorFunction&, const EffectiveScalarFunction&)) &join);
    def("join", (EffectiveVectorFunction(*)(const EffectiveScalarFunction&, const EffectiveVectorFunction&)) &join);
    def("join", (EffectiveVectorFunction(*)(const EffectiveVectorFunction&, const EffectiveVectorFunction&)) &join);

    def("compose", (EffectiveScalarFunction(*)(const EffectiveScalarFunction&,const EffectiveVectorFunction&)) &compose);
    def("compose", (EffectiveVectorFunction(*)(const EffectiveVectorFunction&,const EffectiveVectorFunction&)) &compose);

    from_python<EffectiveVectorFunction>();

    class_<ValidatedVectorFunction> interval_vector_function_class("ValidatedVectorFunction", init<ValidatedVectorFunction>());
    interval_vector_function_class.def(init<EffectiveVectorFunction>());
    interval_vector_function_class.def(init<uint,uint>());
    interval_vector_function_class.def("result_size", &ValidatedVectorFunction::result_size);
    interval_vector_function_class.def("argument_size", &ValidatedVectorFunction::argument_size);
    interval_vector_function_class.def("__getitem__", &ValidatedVectorFunction::get);
    interval_vector_function_class.def("__setitem__", &ValidatedVectorFunction::set);
    interval_vector_function_class.def("__str__", &__cstr__<ValidatedVectorFunction>);
    interval_vector_function_class.def("__repr__", &__crepr__<ValidatedVectorFunction>);

    implicitly_convertible<EffectiveVectorFunction,ValidatedVectorFunction>();
}


void export_scalar_python_function()
{
    class_<ScalarPythonFunction, bases< EffectiveScalarFunctionInterface > > scalar_python_function_class("ScalarUserFunction", init<object>());
    scalar_python_function_class.def(init<uint,object>());
}

void export_vector_python_function()
{
    class_<VectorPythonFunction, bases< EffectiveVectorFunctionInterface > > vector_python_function_class("VectorUserFunction", init<object>());
    vector_python_function_class.def(init<uint,uint,object>());
}


void function_submodule() {
    to_python< Array<std::string> >();
    from_python< Array<std::string> >();

    export_multi_index();

    export_scalar_function();
    export_vector_function();

    //export_scalar_python_function();
    //export_vector_python_function();
}

