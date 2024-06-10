#ifndef ARIADNE_POLYNOMIAL_APPROXIMATION_INTERFACE_HPP
#define ARIADNE_POLYNOMIAL_APPROXIMATION_INTERFACE_HPP

#include "utility/cache.hpp"
#include "utility/factorials.hpp"

namespace Ariadne {
template<typename T>
class IPolynomialApproximation {
  protected:
    using RND = T::RoundingModeType;
    using PR = T::PrecisionType;

  public:
    ~IPolynomialApproximation() {}
    virtual Bounds<T> evaluate(const Bounds<T> &x, int subIntervals = 1) const = 0;
    virtual Bounds<T> evaluateRaw(const Bounds<T> &x) const {
        return evaluate(x);
    }
    virtual Bounds<T> evaluateDerivative(const Bounds<T> &x) const = 0;

    virtual DegreeType degree() const = 0;
    virtual PR precision() const = 0;

    Bounds<T> operator()(const Bounds<T> &x) const {
        return evaluate(x);
    }
};

} // namespace Ariadne

#endif
