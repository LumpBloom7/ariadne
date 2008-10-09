/***************************************************************************
 *            standard_integrator.code.h
 *
 *  Copyright  2006-7  Alberto Casagrande, Pieter Collins
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
 
//#define DEBUG

#include <iosfwd>
#include <string>
#include <sstream>
#include <algorithm>
#include <typeinfo>

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include "standard_flower.h"

#include "base/array.h"
#include "base/exceptions.h"
#include "numeric/interval.h"

#include "linear_algebra/vector.h"
#include "linear_algebra/covector.h"
#include "linear_algebra/matrix.h"
#include "linear_algebra/matrix_function.h"

#include "differentiation/affine_variable.h"
#include "differentiation/differential_vector.h"

#include "function/affine_model.h"
#include "function/taylor_model.h"

#include "geometry/box.h"

#include "system/vector_field.h"

#include "evaluation/bounder_interface.h"

#include "output/logging.h"

#include "differentiation/power_series.code.h"
#include "differentiation/differential.code.h"
#include "differentiation/differential_vector.code.h"

namespace {

using namespace Ariadne;

template<class R>
Matrix<R>
symmetrize(const Vector< Interval<R> >& iv)
{
  Matrix<R> A(iv.size(),iv.size()+1);
  for(size_type i=0; i!=A.number_of_rows(); ++i) {
    A(i,i)=radius(iv(i));
    A(i,iv.size())=midpoint(iv(i));
  }
  return A;
}

} // namespace



namespace Ariadne { 





template<class R>
Point< Interval<R> >
StandardFlower<R>::flow_step(const VectorField<R>& vector_field, 
                                         const Point<I>& initial_point, 
                                         const Rational& step_size, 
                                         const Box<R>& bounding_box) const
{
  typedef Interval<R> I;
  
  // Use second order formula \f$ \Phi(t,p) = p + tf(p) + t^2/2 Df(B)f(B) \f$
  const VectorField<R>& vf=vector_field;
  const Point<I>& p=initial_point;
  Point<I> b=bounding_box;
  I h=step_size;
  
  return p + h * ( vf(p) + (h/2) * ( vf.jacobian(b) * vf(b) ) );
}

template<class R>
Matrix< Interval<R> >
StandardFlower<R>::flow_step_jacobian(const VectorField<R>& vector_field, 
                                                  const Point<I>& initial_point, 
                                                  const Rational& step_size, 
                                                  const Box<R>& bounding_box) const
{
  // Use first order formula \f$ D\Phi(t,p) = I + t Df(B) W \f$ where W is a bound for D\Phi([0,h],p)
  // Use ||W-I|| < e^{Lh}-1, where L is the  norm of Df
  typedef Interval<R> I;
  const dimension_type d=vector_field.dimension();
  const VectorField<R>& vf=vector_field;
  //const Point<I>& pr=initial_point;
  Point<I> b=bounding_box;
  I h=step_size;

  Matrix<I> Id = Matrix<I>::identity(d);

  Matrix<I> Df = vf.jacobian(b);
  R l = norm(Df).upper();
  I e = sub_up(exp_up(mul_up(h.upper(),l)),R(1))*I(-1,1);

  Matrix<I> W = Matrix<I>::identity(d)+e*Matrix<I>::one(d,d);

  // Perform a couple of steps
  W=Id + I(0,h.upper()) * (Df * W);
  W=Id + h * (Df * W);
  return W;
}








template<class X>
class PowerSeriesAffineVariable
  : public PowerSeries< AffineVariable<X> >
{
 public:
  PowerSeriesAffineVariable() : PowerSeries< AffineVariable<X> >() { }
  PowerSeriesAffineVariable(const PowerSeries< AffineVariable<X> >& x)
    : PowerSeries< AffineVariable<X> >(x) { }
  static PowerSeriesAffineVariable<X> constant(uint n, uint d, X v) {
    PowerSeries< AffineVariable<X> > x(d);
    x[0]=AffineVariable<X>::constant(n,v);
    for(uint j=1; j<=d; ++j) {
      x[j]=AffineVariable<X>::constant(n,0.0); 
    }
    return x;
  }
  static PowerSeriesAffineVariable<X> constant_variable(uint n, uint d, X v, uint i) {
    PowerSeries< AffineVariable<X> > x(d);
    x[0]=AffineVariable<X>::variable(n,v,i);
    for(uint j=1; j<=d; ++j) {
      x[j]=AffineVariable<X>::constant(n,0.0); 
    }
    return x;
  }
  static PowerSeriesAffineVariable<X> variable(uint n, uint d, X v, uint i) {
    ARIADNE_ASSERT(d>=1);
    PowerSeries< AffineVariable<X> > x(d);
    x[0]=AffineVariable<X>::variable(n,v,i);
    x[1]=AffineVariable<X>::variable(n,1.0,i);
    for(uint j=2; j<=d; ++j) {
      x[j]=AffineVariable<X>::constant(n,0.0); 
    }
    return x;
  }
};


template<class X>
class PowerSeriesDifferential
  : public PowerSeries< Differential<X> >
{
 public:
  PowerSeriesDifferential() : PowerSeries< Differential<X> >() { }
  PowerSeriesDifferential(const PowerSeries< Differential<X> >& x)
    : PowerSeries< Differential<X> >(x) { }
  static PowerSeriesDifferential<X> constant_variable(uint n, uint ot, uint ox, X v, uint i) {
    PowerSeries< Differential<X> > x(ot);
    x[0]=Differential<X>::variable(n,ox,v,i);
    for(uint j=1; j<=ot; ++j) {
      x[j]=Differential<X>::constant(n,ox,0.0); 
    }
    return x;
  }
  static PowerSeriesDifferential<X> variable(uint n, uint ot, uint ox, X v, uint i) {
    PowerSeries< Differential<X> > x(ot);
    x[0]=Differential<X>::variable(n,ox,v,i);
    x[1]=Differential<X>::variable(n,ox,1.0,i);
    for(uint j=2; j<=ot; ++j) {
      x[j]=Differential<X>::constant(n,ox,0.0); 
    }
    return x;
  }
};

template<class R>
AffineVariable<R> midpoint(const AffineVariable< Interval<R> >& iav)
{
  R v=midpoint(iav.value());
  Covector<R> cv=midpoint(iav.derivative());
  return AffineVariable<R>(v,cv);
}



template<class X> 
array< PowerSeries< AffineVariable<X> > >
integrate(const DifferentialVector<X>& vf, const Point<X> x)
{
  ARIADNE_ASSERT(vf.result_size()==vf.argument_size());
  ARIADNE_ASSERT(vf.argument_size()==x.dimension());
  dimension_type n=x.dimension();
  smoothness_type d=vf.degree();
  array< PowerSeries< AffineVariable<X> > > y(n);
  array< PowerSeries< AffineVariable<X> > > yp(n);
  for(size_type i=0; i!=n; ++i) {
    y[i]=PowerSeries< Differential<X> >(0);
    y[i][0]=AffineVariable<X>::variable(n,x[i],i);
  }
  for(uint j=0; j<d; ++j) {
    yp=evaluate(vf,y);
    //cout << "j="<<j<<"\n y="<<y<<"\n yp="<<yp<<endl;
    for(uint i=0; i!=n; ++i) {  
      y[i]=antiderivative(yp[i],y[i][0]);
    }
  } 
  return y;
}


template<class X> 
array< PowerSeries< Differential<X> > >
integrate(const DifferentialVector<X>& vf, const Point<X> x, smoothness_type ox)
{
  ARIADNE_ASSERT(vf.result_size()==vf.argument_size());
  ARIADNE_ASSERT(vf.argument_size()==x.dimension());
  dimension_type n=x.dimension();
  smoothness_type ot=vf.degree();
  array< PowerSeries< Differential<X> > > y(n);
  array< PowerSeries< Differential<X> > > yp(n);
  for(size_type i=0; i!=n; ++i) {
    y[i]=PowerSeries< Differential<X> >(0);
    y[i][0]=Differential<X>::variable(n,ox,x[i],i);
  }
  for(uint j=0; j<ot; ++j) {
    yp=evaluate(vf,y);
    //cout << "j="<<j<<"\n y="<<y<<"\n yp="<<yp<<endl;
    for(uint i=0; i!=n; ++i) {  
      y[i]=antiderivative(yp[i],y[i][0]);
    }
  } 
  return y;
}
















template<class R>
AffineModel<R>
StandardFlower<R>::affine_flow_model(const VectorField<R>& vector_field, 
                                     const Point<R>& initial_point, 
                                     const Box<R>& initial_domain, 
                                     const Rational& step_size, 
                                     const Box<R>& bounding_box) const
{
  // Convert from Differential flow model
  
  typedef Interval<R> I;

  TaylorModel<R> taylor_flow_model=this->taylor_flow_model(vector_field,initial_point,initial_domain,step_size,bounding_box);
  return AffineModel<R>(taylor_flow_model.domain(),
                        taylor_flow_model.centre(),
                        taylor_flow_model.evaluate(taylor_flow_model.centre()),
                        taylor_flow_model.jacobian(taylor_flow_model.domain()));

  uint to=this->temporal_order();
  dimension_type n=initial_point.dimension();
  I h=step_size;
  
  // Make dvf contain the vector field derivatives at the centre of the initial set,
  // except for the highest-order term, which contains the derivatives over the entire set.
  DifferentialVector<I> dvf=vector_field.derivative(Point<I>(bounding_box),to);
  DifferentialVector<I> cvf=vector_field.derivative(initial_point,to-1);
  for(uint i=0; i!=n; ++i) {
    dvf[i].assign(cvf[i]);
  }
 
  // Set up array of flow derivative values
  // Each component is a constant in time and a variable in space.
  array< PowerSeriesAffineVariable<I> > y(n);
  for(size_type i=0; i!=n; ++i) {
    y[i]=PowerSeriesAffineVariable<I>::constant_variable(n,0,initial_point[i],i);
  }

  // Compute the Differential series of the state and first variation
  array< PowerSeriesAffineVariable<I> > yp(n);
  for(uint j=0; j<to; ++j) {
    yp=evaluate(dvf,y);
    for(uint i=0; i!=n; ++i) {  
      y[i]=antiderivative(yp[i],y[i][0]);
    }
  } 

  //for(uint j=0; j<=to; ++j) { cout << "y["<<j<<"]=\n"; for(uint i=0; i!=n; ++i) { cout << " " << midpoint(y[i][j]) << endl; } }

  // Compute the state and first variation at the final time
  Vector< AffineVariable<I> > r(n);
  for(uint i=0; i!=n; ++i) {
    r[i]=y[i][0];
  }
  I c=1;
  for(uint j=1; j<=to; ++j) {
    c*=h; c/=j;
    for(uint i=0; i!=n; ++i) {
      r[i] += y[i][j]*c;
    }
  }

  return AffineModel<R>(bounding_box.position_vectors(),initial_point.position_vector(),r);
}



template<class R>
TaylorModel<R>
StandardFlower<R>::taylor_flow_model(const VectorField<R>& vector_field, 
                                                 const Point<R>& initial_point, 
                                                 const Box<R>& initial_domain, 
                                                 const Rational& step_size, 
                                                 const Box<R>& bounding_box) const
{
  uint verbosity=0;
  ARIADNE_LOG(6,"taylor_flow_model(...)\n");
  
  dimension_type n=initial_point.dimension();
  smoothness_type ot=this->temporal_order();
  smoothness_type ox=this->spacial_order();
  Interval<R> h=step_size;

  DifferentialVector<I> vfc=vector_field.derivative(initial_point,ot);
  DifferentialVector<I> vfb=vector_field.derivative(bounding_box,ot);

  const array<R>& x=initial_point.data();
  array< PowerSeries< Differential<I> > > y(n);
  for(uint i=0; i!=n; ++i) {
    // y[i][0]=Differential<I>::variable(n,ox,x[i],i);
    y[i][0]=Differential<I>::variable(n,ox,0.0,i);
  }
  ARIADNE_LOG(7,"y="<<y<<"\n");
  array< PowerSeries< Differential<I> > > yp(n);
  for(uint j=0; j<=ot; ++j) {
    yp=evaluate(vfc,y);
    for(uint i=0; i!=n; ++i) {  
      y[i]=antiderivative(yp[i],y[i][0]);
    }
    ARIADNE_LOG(7,"yp="<<yp<<"\ny="<<y<<"\n");
  }

  for(uint i=0; i!=n; ++i) {  
    y[i][0].value()=x[i];
  }
  ARIADNE_LOG(7,"\ny="<<y<<"\n\n");

  DifferentialVector<I> phi(n,n,ox);
  for(uint j=0; j<=ot; ++j) {
    for(uint i=0; i!=n; ++i) {
      phi[i]+=y[i][j]*pow(h,j);
    }
    ARIADNE_LOG(7,"phi="<<phi<<"\n");
  }

  //FIXME: Put rigorous error bounds in flow model
  return TaylorModel<R>(initial_domain.position_vectors(),initial_point.position_vector(),phi,phi);
}




} // namespace Ariadne