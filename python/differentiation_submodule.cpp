/***************************************************************************
 *            differentiation_submodule.cpp
 *
 *  Copyright 2008--17  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "boost_python.hpp"
#include "utilities.hpp"

#include "utility/typedefs.hpp"
#include "utility/array.hpp"
#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "algebra/covector.hpp"
#include "algebra/differential.hpp"
#include "algebra/univariate_differential.hpp"
#include "algebra/fixed_differential.hpp"
#include "algebra/fixed_univariate_differential.hpp"
#include "algebra/expansion.inl.hpp"

using namespace boost::python;
using namespace Ariadne;

template<class X> Void read_array(Array<X>&, const boost::python::object& obj) { }
inline Nat compute_polynomial_data_size(Nat rs, Nat as, Nat d) { return rs*Ariadne::bin(d+as,as); }

namespace Ariadne {

// FIXME: Ensure all valid arithmetic and comparisons are defined!
inline auto operator==(FloatDPBounds x, Int n) -> decltype(x==FloatDPValue(n)) { return x==FloatDPValue(n); }
inline auto operator!=(FloatDPBounds x, Int n) -> decltype(x!=FloatDPValue(n)) { return x!=FloatDPValue(n); }
inline auto operator> (FloatDPBounds x, Int n) -> decltype(x> FloatDPValue(n)) { return x> FloatDPValue(n); }
inline auto operator*=(FloatDPApproximation x, Int n) -> decltype(x*=FloatDPApproximation(n)) { return x*=FloatDPApproximation(n); }

template<class X>
struct to_python_dict< Ariadne::Expansion<MultiIndex,X>  > {
    to_python_dict() { boost::python::to_python_converter< Ariadne::Expansion<MultiIndex,X>, to_python_dict< Ariadne::Expansion<MultiIndex,X> > >(); }
    static PyObject* convert(const Ariadne::Expansion<MultiIndex,X>& e) {
        Nat n=e.argument_size();
        boost::python::dict res;
        boost::python::list lst;
        for(Nat i=0; i!=n; ++i) { lst.append(0); }
        Ariadne::MultiIndex a;
        X c;
        for(typename Expansion<MultiIndex,X>::ConstIterator iter=e.begin(); iter!=e.end(); ++iter) {
            a=iter->index();
            c=iter->coefficient();
            for(Nat i=0; i!=a.size(); ++i) { Int ai=a[i]; lst[i]=ai; }
            boost::python::tuple tup(lst);
            //res[tup]=boost::python::object(c);
            res[boost::python::object(a)]=boost::python::object(c);
        }
        return boost::python::incref(boost::python::dict(res).ptr());
    }
    static const PyTypeObject* get_pytype() { return &PyDict_Type; }
};

template<class X>
struct to_python_list< Ariadne::Expansion<MultiIndex,X>  > {
    to_python_list() { boost::python::to_python_converter< Ariadne::Expansion<MultiIndex,X>, to_python_list< Ariadne::Expansion<MultiIndex,X> > >(); }
    static PyObject* convert(const Ariadne::Expansion<MultiIndex,X>& e) {
        Nat n=e.argument_size();
        boost::python::list res;
        boost::python::list alst;
        for(Nat i=0; i!=n; ++i) { alst.append(0); }
        std::cerr<<"Here\n";
        boost::python::list pr; pr.append(0); pr.append(0);
        Ariadne::MultiIndex a;
        X c;
        for(typename Expansion<MultiIndex,X>::ConstIterator iter=e.begin(); iter!=e.end(); ++iter) {
            a=iter->index();
            c=iter->coefficient();
            for(Nat i=0; i!=n; ++i) { Int ai=a[i]; alst[i]=ai; }
            pr[0]=boost::python::tuple(alst);
            pr[1]=boost::python::object(c);
            res.append(boost::python::tuple(pr));
        }
        return boost::python::incref(boost::python::list(res).ptr());
    }
    static const PyTypeObject* get_pytype() { return &PyList_Type; }
};

/*
   static PyObject* convert(const Tuple<T1,T2,T3,T4>& tup) {
        boost::python::list lst;
        lst.append(boost::python::object(tup.first));
        lst.append(boost::python::object(tup.second));
        lst.append(boost::python::object(tup.third));
        lst.append(boost::python::object(tup.fourth));
        boost::python::tuple result(lst);
        return boost::python::incref(boost::python::tuple(result).ptr());
    static PyObject* convert(const std::map<K,V>& map) {
        boost::python::dict result;
        for(typename std::map<K,V>::ConstIterator iter=map.begin(); iter!=map.end(); ++iter) {
            result[boost::python::object(iter->first)]=boost::python::object(iter->second);
        }
        return boost::python::incref(boost::python::dict(result).ptr());
*/
}


template<class DIFF>
DIFF*
make_dense_differential(const Nat& as, const Nat& d, const boost::python::object& obj)
{
    typedef typename DIFF::ValueType X;
    DIFF* result=new DIFF(as,d);
    Array<X> coefficient;
    read_array(coefficient,obj);
    std::cerr<<"polynomial_data_size("<<as<<","<<d<<")="<<compute_polynomial_data_size(1u,as,d)<<"\n";
    std::cerr<<"coefficient="<<coefficient<<"\n";
    assert(coefficient.size()==compute_polynomial_data_size(1u,as,d));
    MultiIndex i(as);
    const X* ptr=coefficient.begin();
    while(i.degree()<=d) {
        //result[i]=*ptr; ++i; ++ptr;
        result->expansion().append(i,*ptr); ++i; ++ptr;
    }
    return result;
}

template<class DIFF>
DIFF*
make_sparse_differential(const boost::python::object& obj,const Nat& d)
{
    typedef typename DIFF::ValueType X;
    Expansion<MultiIndex,X> expansion = boost::python::extract< Expansion<MultiIndex,X> >(obj);
    DIFF* result=new DIFF(expansion,d);
    return result;
}


template<class DIFF>
boost::python::list
make_differential_variables(const Nat& d, const Vector<typename DIFF::NumericType>& x)
{
    boost::python::list result;
    for(Nat i=0; i!=x.size(); ++i) {
        result.append(DIFF::variable(x.size(),d,x[i],i));
    }
    return result;
}


template<class DIFF>
Vector<DIFF>*
make_differential_vector(const Nat& rs, const Nat& as, const Nat& d, const boost::python::object& obj)
{
    typedef typename DIFF::ValueType X;
    boost::python::list lst=boost::python::extract<boost::python::list>(obj);
    Array<X> coefficient;
    read_array(coefficient,obj);
    ARIADNE_ASSERT(coefficient.size()==compute_polynomial_data_size(rs,as,d));
    Vector<DIFF>* result=new Vector<DIFF>(rs,DIFF(as,d));
    for(SizeType i=0; i!=rs; ++i) {
        DIFF* ri = make_sparse_differential<DIFF>(lst[i],d);
        (*result)[i]=*ri;
        delete ri;
    }
    return result;
}

template<class PR> Differential<FloatBounds<PR>> make_differential_variable(SizeType n, DegreeType d, Real r, SizeType i, PR pr) {
    return Differential<FloatBounds<PR>>::variable(n,d,FloatBounds<PR>(r,pr),i);
}
template<class PR> Differential<FloatBounds<PR>> make_differential_variable(SizeType n, DegreeType d, ValidatedReal r, SizeType i, PR pr) {
    return Differential<FloatBounds<PR>>::variable(n,d,FloatBounds<PR>(r,pr),i);
}
template<class X, class Y, class PR> Vector<Differential<X>> make_variables(DegreeType d, Vector<Y> vy, PR pr) {
    return Differential<X>::variables(d,Vector<X>(vy,pr));
}


template<class C, class I, class X, EnableIf<IsSame<I,Int>> =dummy> inline
X get_item(const C& c, const I& i) {
    return c[static_cast<Nat>(i)];
}

template<class C, class I, class X, EnableIf<Not<IsSame<I,Int>>> =dummy> inline
X get_item(const C& c, const I& i) {
    return c[i];
}

template<class C, class I, class J, class X> inline
X matrix_get_item(const C& c, const I& i, const J& j) { return c[static_cast<SizeType>(i)][j]; }

template<class C, class I, class X, EnableIf<IsSame<I,Int>> =dummy> inline
Void set_item(C& c, const I& i, const X& x) { c[static_cast<SizeType>(i)]=x; }

template<class C, class I, class X, EnableIf<Not<IsSame<I,Int>>> =dummy> inline
Void set_item(C& c, const I& i, const X& x) { c[i]=x; }

template<class C, class I, class J, class X> inline
Void matrix_set_item(C& c, const I& i, const J& j, const X& x) { c[i][j]=x; }


namespace Ariadne {

template<class X> OutputStream& operator<<(OutputStream& os, const PythonRepresentation< Expansion<MultiIndex,X> >& repr);

template<class X> OutputStream& operator<<(OutputStream& os, const PythonRepresentation< Differential<X> >& repr) {
    const Differential<X>& diff=repr.reference();
    os << python_name<X>("Differential").c_str() << "(" << python_representation(diff.expansion()) << "," << diff.degree() << ")";
    //os << python_name<X>("Differential").c_str() << "(" << diff.argument_size() << "," << diff.degree() << "," << python_representation(diff.expansion()) << ")";
    return os;
}

}


template<class D> using ValueType = decltype(declval<D>().value());
template<class D> using GradientType = decltype(declval<D>().gradient());
template<class D> using HessianType = decltype(declval<D>().hessian());

template<class DIFF>
Void export_differential(const String& name)
{
    typedef typename DIFF::ValueType X;
    typedef typename X::GenericType Y;
    typedef DIFF D;

    class_<D> differential_class(name.c_str(), init<D>() );
    differential_class.def("__init__", make_constructor(&make_sparse_differential<D>) );
    differential_class.def( init< SizeType, DegreeType >());
    differential_class.def( init< Expansion<MultiIndex,X>, DegreeType >());
    differential_class.def("__getitem__", &get_item<D,MultiIndex,X>);
    differential_class.def("__setitem__",&set_item<D,MultiIndex,X>);
    differential_class.def(-self);
    differential_class.def(self+self);
    differential_class.def(self-self);
    differential_class.def(self*self);
    differential_class.def(self/self);
    differential_class.def(self+X());
    differential_class.def(self-X());
    differential_class.def(self*X());
    differential_class.def(self/X());
    differential_class.def(X()+self);
    differential_class.def(X()-self);
    differential_class.def(X()*self);
    differential_class.def(self+Y());
    differential_class.def(self-Y());
    differential_class.def(self*Y());
    differential_class.def(self/Y());
    differential_class.def(Y()+self);
    differential_class.def(Y()-self);
    differential_class.def(Y()*self);
    differential_class.def(self+=self);
    differential_class.def(self-=self);
    differential_class.def(self+=X());
    differential_class.def(self-=X());
    differential_class.def(self*=X());
    differential_class.def(self/=X());
    differential_class.def(self_ns::str(self));
    differential_class.def("__repr__", &__repr__<D>);

    differential_class.def("value", (ValueType<D>(D::*)()const)&D::value, return_value_policy<copy_const_reference>());
    differential_class.def("gradient", (GradientType<D>(D::*)()const)&D::gradient);
    differential_class.def("hessian", (HessianType<D>(D::*)()const)&D::hessian);
    differential_class.def("expansion", (Expansion<MultiIndex,X>const&(D::*)()const)&D::expansion, return_value_policy<copy_const_reference>());

    differential_class.def("constant",(D(*)(SizeType, DegreeType, const X&))&D::constant);
    differential_class.def("variable",(D(*)(SizeType, DegreeType, const X&, SizeType))&D::variable);
    differential_class.def("constants",(Vector<D>(*)(SizeType, DegreeType, const Vector<X>&))&D::constants);
    differential_class.def("variables",(Vector<D>(*)(DegreeType, const Vector<X>&))&D::variables);


    differential_class.staticmethod("constant");
    differential_class.staticmethod("variable");
    differential_class.staticmethod("variables");

    def("derivative", (D(*)(const D&, SizeType))&D::_derivative);
    def("antiderivative", (D(*)(const D&, SizeType))&D::_antiderivative);

    def("neg",&_neg_<D>);
    def("sqr",&_sqr_<D>);
    def("rec",&_rec_<D>);
    def("pow",&_pow_<D,Int>);
    def("sqrt",&_sqrt_<D>);
    def("exp",&_exp_<D>);
    def("log",&_log_<D>);
    def("sin",&_sin_<D>);
    def("cos",&_cos_<D>);
    def("tan",&_tan_<D>);
    def("atan",&_atan_<D>);
}

template<class DIFF>
Void
export_differential_vector(const String& name)
{
    typedef typename DIFF::ValueType X;
    typedef Vector<X> V;
    typedef DIFF D;
    typedef Vector<D> DV;

    class_<DV> differential_vector_class(name.c_str(), init<DV>());
    differential_vector_class.def("__init__", make_constructor(&make_differential_vector<D>) );
    differential_vector_class.def( init< Nat, Nat, Nat >());
    differential_vector_class.def("__getitem__", &matrix_get_item<DV,Int,MultiIndex,X>);
    differential_vector_class.def("__getitem__", &get_item<DV,Int,D>);
    differential_vector_class.def("__setitem__",&set_item<DV,Int,X>);
    differential_vector_class.def("__setitem__",&set_item<DV,Int,D>);
    differential_vector_class.def("__neg__",&__neg__<DV,DV>);
    differential_vector_class.def("__add__",&__add__<DV,DV,DV>);
    differential_vector_class.def("__sub__",&__sub__<DV,DV,DV>);
    differential_vector_class.def("__add__",&__add__<DV,DV,V>);
    differential_vector_class.def("__sub__",&__sub__<DV,DV,V>);
    //differential_vector_class.def("__mul__",&__mul__<DV,DV,D>);
    //differential_vector_class.def("__div__",&__div__<DV,DV,D>);
    differential_vector_class.def("__rmul__",&__rmul__<DV,DV,X>);
    differential_vector_class.def("__mul__",&__mul__<DV,DV,X>);
    differential_vector_class.def("__div__",&__div__<DV,DV,X>);
    differential_vector_class.def("value", &DV::value);
    differential_vector_class.def("jacobian", &DV::jacobian);
    differential_vector_class.def(self_ns::str(self));
    differential_vector_class.def(self_ns::repr(self));

    def("compose",(D(*)(const D&,const DV&))&D::_compose);
    def("compose",(DV(*)(const DV&,const DV&))&DV::_compose);

    def("solve",(DV(*)(const DV&,const V&))&DV::_solve);
    def("flow",(DV(*)(const DV&,const V&))&DV::_flow);

    //def("lie_derivative", (DV(*)(const DV&,const DV&))&lie_derivative);
}

template Void export_differential< Differential<FloatDPApproximation> >(const String&);
template Void export_differential< Differential<FloatDPBounds> >(const String&);

template Void export_differential_vector< Differential<FloatDPApproximation> >(const String&);
template Void export_differential_vector< Differential<FloatDPBounds> >(const String&);

Void differentiation_submodule()
{
    to_python_dict < Expansion<MultiIndex,FloatDPApproximation> >();
    to_python_dict < Expansion<MultiIndex,FloatDPBounds> >();

    export_differential< Differential<FloatDPApproximation> >(python_name<FloatDPApproximation>("Differential"));
    export_differential< Differential<FloatDPBounds> >(python_name<FloatDPBounds>("Differential"));
    export_differential_vector< Differential<FloatDPApproximation> >(python_name<FloatDPApproximation>("DifferentialVector"));
    export_differential_vector< Differential<FloatDPBounds> >(python_name<FloatDPBounds>("DifferentialVector"));

    export_differential< Differential<FloatMPApproximation> >(python_name<FloatMPApproximation>("Differential"));
    export_differential< Differential<FloatMPBounds> >(python_name<FloatMPBounds>("Differential"));
    export_differential_vector< Differential<FloatMPApproximation> >(python_name<FloatMPApproximation>("DifferentialVector"));
    export_differential_vector< Differential<FloatMPBounds> >(python_name<FloatMPBounds>("DifferentialVector"));

    def("differential_variables", (Vector<Differential<FloatDPBounds>>(*)(DegreeType, Vector<Real>, DoublePrecision)) &make_variables<FloatDPBounds>);
    def("differential_variables", (Vector<Differential<FloatMPBounds>>(*)(DegreeType, Vector<Real>, MultiplePrecision)) &make_variables<FloatMPBounds>);


}

