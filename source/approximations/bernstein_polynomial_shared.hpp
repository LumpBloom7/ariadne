#ifndef ARIADNE_BERNSTEIN_POLYNOMIAL_SHARED_HPP
#define ARIADNE_BERNSTEIN_POLYNOMIAL_SHARED_HPP

#include "utility/cache.hpp"
#include "utility/factorials.hpp"

namespace Ariadne {
class IBernsteinPolynomialBase {
  protected:
    static SimpleCache<Integer, Nat64, Nat64> binomialCoefficients;
};

template<typename T>
class IBernsteinPolynomial {
  protected:
    using RND = T::RoundingModeType;
    using PR = T::PrecisionType;

  public:
    ~IBernsteinPolynomial() {}
    virtual Bounds<T> evaluate(const Bounds<T> &x) const = 0;
    virtual Bounds<T> evaluateDerivative(const Bounds<T> &x) const = 0;

    virtual DegreeType degree() const = 0;
    virtual PR precision() const = 0;

    Bounds<T> operator()(const Bounds<T> &x) const {
        return evaluate(x);
    }
};

} // namespace Ariadne

#endif
