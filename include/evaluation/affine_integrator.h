/***************************************************************************
 *            affine_integrator.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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
 
/*! \file affine_integrator.h
 *  \brief Methods for integrating points and sets under a vector field.
 */

#ifndef ARIADNE_AFFINE_INTEGRATOR_H
#define ARIADNE_AFFINE_INTEGRATOR_H

#include "../linear_algebra/declarations.h"
#include "../geometry/declarations.h"
#include "../system/declarations.h"
#include "../evaluation/integrator_interface.h"

namespace Ariadne {
  namespace Evaluation {
   
    /*! \brief Compute an over-approximation to \f$\sum_{n=0}^{\infty} \frac{x^n}{(k+n)!}\f$. */
    template<class R> R gexp_up(const R& x, uint k);
  
    /*! \brief Compute an interval approximation to \f$\sum_{n=0}^{\infty} \frac{(tA)^n}{(k+n)!}b\f$. */
    template<class R> 
    LinearAlgebra::Vector< Numeric::Interval<R> > 
    gexp(const LinearAlgebra::Matrix<R>& A, const LinearAlgebra::Vector<R>& b, const Numeric::Interval<R>& t, const uint& k);
      
    /*! \brief Compute an interval approximation to \f$\sum_{n=0}^{\infty} \frac{(tA)^n}{(k+n)!}\f$. */
    template<class R> 
    LinearAlgebra::Matrix< Numeric::Interval<R> > 
    gexp(const LinearAlgebra::Matrix<R>& Ab, const Numeric::Interval<R>& t, const uint& k);
      
    
    /*!\ingroup Integrate
     * \brief An integrator based on using the exponential formula to integrate an affine vector field. 
     *  
     * The \f$C^1\f$-Affine algorithm is a Taylor method.
     */
    template<class R>
    class AffineIntegrator
      : public DifferentiableIntegratorInterface<R> 
    {
      typedef Numeric::Interval<R> I;
      typedef Geometry::Zonotope<I> BS;
     public:
      /*! \brief Constructor. */
      AffineIntegrator();

      /*! \brief Cloning operator. */
      virtual AffineIntegrator<R>* clone() const;

     public:

      /*! \brief Integrate a basic set for within a bounding set. */
      virtual Geometry::Point<I> flow_step(const System::VectorFieldInterface<R>& vf,
                                           const Geometry::Point<I>& p,
                                           const Numeric::Interval<R>& t,
                                           const Geometry::Rectangle<R>& bb) const;
     
      /*! \brief Integrate a basic set for within a bounding set. */
      virtual LinearAlgebra::Matrix<I> flow_step_jacobian(const System::VectorFieldInterface<R>& vf,
                                                          const Geometry::Point<I>& p,
                                                          const Numeric::Interval<R>& r,
                                                          const Geometry::Rectangle<R>& bb) const;
     

      virtual Geometry::Rectangle<R> integration_step(const System::VectorFieldInterface<R>& vector_field,
                                                     const Geometry::Rectangle<R>& initial_set,
                                                     const Numeric::Interval<R>& step_size,
                                                     const Geometry::Rectangle<R>& bounding_set) const;
      
      virtual Geometry::Rectangle<R> reachability_step(const System::VectorFieldInterface<R>& affine_vector_field,
                                                      const Geometry::Rectangle<R>& initial_set,
                                                      const Numeric::Interval<R>& step_size,
                                                      const Geometry::Rectangle<R>& bounding_set) const;

      virtual Geometry::Zonotope<I,R> integration_step(const System::VectorFieldInterface<R>& vector_field,
                                                     const Geometry::Zonotope<I,R>& initial_set,
                                                     const Numeric::Interval<R>& step_size,
                                                     const Geometry::Rectangle<R>& bounding_set) const;
      
      virtual Geometry::Zonotope<I,R> reachability_step(const System::VectorFieldInterface<R>& affine_vector_field,
                                                        const Geometry::Zonotope<I,R>& initial_set,
                                                      const Numeric::Interval<R>& step_size,
                                                      const Geometry::Rectangle<R>& bounding_set) const;

      /*! \brief A \f$C^\infty\f$ algorithm for integrating forward a zonotope. */
      virtual Geometry::Zonotope<I,I> integration_step(const System::VectorFieldInterface<R>& vector_field,
                                                     const Geometry::Zonotope<I>& initial_set,
                                                     const Numeric::Interval<R>& step_size,
                                                     const Geometry::Rectangle<R>& bounding_set) const;
      
      /*! \brief A \f$C^\infty\f$ algorithm for integrating forward a zonotope for a time up to time \a step_size. */
      virtual Geometry::Zonotope<I,I> reachability_step(const System::VectorFieldInterface<R>& affine_vector_field,
                                                      const Geometry::Zonotope<I>& initial_set,
                                                      const Numeric::Interval<R>& step_size,
                                                      const Geometry::Rectangle<R>& bounding_set) const;
     public:
      /*! \brief Integrate \a initial point for time \a step_size. */
      Geometry::Point<I> flow_step(const System::AffineVectorField<R>& affine_vector_field,
                                              const Geometry::Point<I>& initial_point,
                                              const Numeric::Interval<R>& step_size) const;

      /*! \brief Comput the spacial Jacobian derivative at \a initial point for time \a step_size. */
      LinearAlgebra::Matrix<I> flow_step_jacobian(const System::AffineVectorField<R>& vector_field,
                                                   const Geometry::Point<I>& affine_initial_point,
                                                   const Numeric::Interval<R>& step_size) const;

      /*! \brief A \f$C^\infty\f$ algorithm for integrating forward a zonotope for a time up to time \a step_size. */
      Geometry::Zonotope<I,R> integration_step(const System::AffineVectorField<R>& affine_vector_field,
                                               const Geometry::Zonotope<I,R>& initial_set,
                                               const Numeric::Interval<R>& step_size) const;

      /*! \brief A \f$C^\infty\f$ algorithm for integrating forward a zonotope for a time up to time \a step_size. */
     Geometry::Zonotope<I,R> reachability_step(const System::AffineVectorField<R>& affine_vector_field,
                                               const Geometry::Zonotope<I,R>& initial_set,
                                               const Numeric::Interval<R>& step_size) const;

      /*! \brief A \f$C^\infty\f$ algorithm for integrating forward a zonotope for a time up to time \a step_size. */
      Geometry::Zonotope<I> integration_step(const System::AffineVectorField<R>& affine_vector_field,
                                              const Geometry::Zonotope<I>& initial_set,
                                              const Numeric::Interval<R>& step_size) const;

      /*! \brief A \f$C^\infty\f$ algorithm for integrating forward a zonotope for a time up to time \a step_size. */
     Geometry::Zonotope<I> reachability_step(const System::AffineVectorField<R>& affine_vector_field,
                                               const Geometry::Zonotope<I>& initial_set,
                                               const Numeric::Interval<R>& step_size) const;


     };

     
  }
}

#endif /* ARIADNE_AFFINE_INTEGRATOR_H */
