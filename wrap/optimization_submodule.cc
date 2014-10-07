/***************************************************************************
 *            optimization_submodule.cc
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

#include "boost_python.h"
#include "utilities.h"

#include "config.h"

#include <boost/python.hpp>

#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "box.h"

#include "function.h"
#include "constraint.h"
#include "linear_programming.h"
#include "nonlinear_programming.h"
#include "constraint_solver.h"

using namespace boost::python;
using namespace Ariadne;


template<class X>
boost::python::tuple
python_compute_basis(const Matrix<X>& A) {
    Array<size_t> p;
    Matrix<X> B;
    make_lpair(p,B)=compute_basis(A);
    boost::python::list l;
    for(size_t i=0; i!=p.size(); ++i) {
        l.append(p[i]);
    }
    return boost::python::make_tuple(l,B);
}

template<class T> T get(const Array<T>& ary, size_t i) { return ary[i]; }
template<class T> void set(Array<T>& ary, size_t i, const T& t) { ary[i]=t; }

template<class T>
void export_internal_array(const char* name)
{
    class_< Array<T> > array_class(name,no_init);
    array_class.def("__len__", &Array<T>::size);
    array_class.def("__getitem__",&get<T>);
    array_class.def(boost::python::self_ns::str(self));
}


void export_variable_type()
{
    typedef Array<Slackness> SlacknessArray;

    enum_<Slackness> variable_enum("Slackness");
    variable_enum.value("BASIS", BASIS);
    variable_enum.value("LOWER", LOWER);
    variable_enum.value("UPPER", UPPER);
}

void export_constraint()
{
    class_<EffectiveConstraint> effective_nonlinear_constraint_class("EffectiveConstraint",init<Real,EffectiveScalarFunction,EffectiveNumberType>());
    effective_nonlinear_constraint_class.def(self_ns::str(self));

    class_<ValidatedConstraint> validated_nonlinear_constraint_class("ValidatedConstraint",init<ValidatedNumberType,ValidatedScalarFunction,ValidatedNumberType>());
    validated_nonlinear_constraint_class.def(init<ValidatedConstraint>());
    validated_nonlinear_constraint_class.def(init<EffectiveConstraint>());
    validated_nonlinear_constraint_class.def("lower_bound", &ValidatedConstraint::lower_bound, return_value_policy<copy_const_reference>());
    validated_nonlinear_constraint_class.def("upper_bound", &ValidatedConstraint::upper_bound, return_value_policy<copy_const_reference>());
    validated_nonlinear_constraint_class.def("function", (const ValidatedScalarFunction&(ValidatedConstraint::*)()const) &ValidatedConstraint::function, return_value_policy<copy_const_reference>());
    validated_nonlinear_constraint_class.def(self_ns::str(self));
}

void export_interior_point_solver()
{
    to_python< Ariadne::Tuple< Vector<Float>, Vector<Float>, Vector<Float> > >();

    class_<InteriorPointSolver> interior_point_solver_class("InteriorPointSolver",init<>());
    interior_point_solver_class.def("minimise", &InteriorPointSolver::minimise);
    interior_point_solver_class.def("feasible", (tribool(InteriorPointSolver::*)(const Vector<Float>&,const Vector<Float>&, const Matrix<Float>&,const Vector<Float>&)const) &InteriorPointSolver::feasible);
    interior_point_solver_class.def("validate_feasibility", &InteriorPointSolver::validate_feasibility);
}


void export_constraint_solver()
{
    class_<ConstraintSolver> constraint_solver_class("ConstraintSolver", init<>());
    constraint_solver_class.def("hull_reduce", (bool(ConstraintSolver::*)(UpperBox&,const ValidatedScalarFunctionInterface&,const ExactInterval&)const) &ConstraintSolver::hull_reduce);
    constraint_solver_class.def("box_reduce", (bool(ConstraintSolver::*)(UpperBox&,const ValidatedScalarFunctionInterface&,const ExactInterval&,uint)const) &ConstraintSolver::box_reduce);
    constraint_solver_class.def("monotone_reduce", (bool(ConstraintSolver::*)(UpperBox&,const ValidatedScalarFunctionInterface&,const ExactInterval&,uint)const) &ConstraintSolver::monotone_reduce);
    constraint_solver_class.def("reduce", (bool(ConstraintSolver::*)(UpperBox&,const List<ValidatedConstraint>&)const) &ConstraintSolver::reduce);
    constraint_solver_class.def("reduce", (bool(ConstraintSolver::*)(UpperBox&,const ValidatedVectorFunction&,const ExactBox&)const) &ConstraintSolver::reduce);
}



template<class X>
void export_simplex_solver()
{
    typedef Array<size_t> SizeArray;

    to_python< std::pair< Array<size_t>, Matrix<X> > >();

    class_< SimplexSolver<X> > simplex_solver_class("SimplexSolver", init<>());
    simplex_solver_class.def("lpstep",(bool(SimplexSolver<X>::*)(const Vector<X>&,const Vector<X>&,const Vector<X>&,const Matrix<X>&,const Vector<X>&,Array<Slackness>& ,SizeArray&,Matrix<X>&,Vector<X>&)const) &SimplexSolver<X>::lpstep);


    simplex_solver_class.def("feasible",(tribool(SimplexSolver<X>::*)(const Vector<X>&,const Vector<X>&,const Matrix<X>&,const Vector<X>&)const) &SimplexSolver<X>::feasible);

    simplex_solver_class.def("verify_feasibility",(tribool(SimplexSolver<X>::*)(const Vector<X>&,const Vector<X>&,const Matrix<X>&,const Vector<X>&,const Array<Slackness>&)const) &SimplexSolver<X>::verify_feasibility);

    simplex_solver_class.def("compute_basis",(std::pair< SizeArray, Matrix<X> >(SimplexSolver<X>::*)(const Matrix<X>&)const) &SimplexSolver<X>::compute_basis);

}


void optimization_submodule() {
    export_variable_type();
    export_constraint();
    export_array<size_t>("SizeArray");
    export_internal_array<Slackness>("SlacknessArray");
    export_simplex_solver<Float>();
#ifdef HAVE_GMPXX_H
    export_simplex_solver<Rational>();
#endif
    export_interior_point_solver();
    export_constraint_solver();
}
