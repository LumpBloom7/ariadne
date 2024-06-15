#ifndef ARIADNE_POLYNOMIAL_APPROXIMATION_INTERFACE_HPP
#define ARIADNE_POLYNOMIAL_APPROXIMATION_INTERFACE_HPP

#include "geometry/interval.hpp"
#include "utility/cache.hpp"
#include "utility/factorials.hpp"

namespace Ariadne {
template<typename T> class IPolynomialApproximation {
  protected:
    using RND = T::RoundingModeType;
    using PR = T::PrecisionType;

  public:
    ~IPolynomialApproximation() {}
    virtual Bounds<T> apply(const Interval<T> &x, Nat subIntervals = 1) const = 0;
    virtual Bounds<T> evaluate(const Bounds<T> &x) const { return apply(toInterval(x)); }
    virtual Bounds<T> evaluateDerivative(const Bounds<T> &x) const = 0;

    virtual DegreeType degree() const = 0;
    virtual PR precision() const = 0;

    Bounds<T> operator()(const Bounds<T> &x) const { return evaluate(x); }
    Bounds<T> operator()(const Interval<T> &x, Nat subintervals = 1) const { return apply(x, subintervals); }

  protected:
    Interval<T> toInterval(const Bounds<T> &bounds) const { return Interval<T>(bounds._l, bounds._u); }
    Bounds<T> fromInterval(const Interval<T> &interval)const  { return cast_singleton(interval); }
};

} // namespace Ariadne

#endif
