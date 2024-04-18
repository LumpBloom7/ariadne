#ifndef ARIADNE_BERNSTEIN_POLYNOMIAL_HPP
#define ARIADNE_BERNSTEIN_POLYNOMIAL_HPP

#include <cmath>
#include <functional>
#include <iostream>

#include "numeric/float_bounds.hpp"
#include "numeric/floatmp.hpp"
#include "numeric/integer.hpp"
#include "utility/cache.hpp"
#include "utility/factorials.hpp"
#include "utility/hash.hpp"
#include "utility/hash_numeric.hpp"
#include "utility/hash_tuple.hpp"
#include "utility/standard.hpp"
#include "numeric/float_approximation.hpp"

namespace Ariadne {

class BernsteinPolynomialBase {
  protected:
    static SimpleCache<Integer, Nat64, Nat64> binomialCoefficients;
};

template<typename T>
    requires HasPrecisionType<T> && ConvertibleTo<T, Bounds<T>>
class BernsteinPolynomial : protected BernsteinPolynomialBase {
  public:
    BernsteinPolynomial(const std::function<T(T)>& function, int degree, size_t precision) : _function{function}, _degree{degree}, _precision{precision} {
        _coefficients = computeCoefficients();
    }

    Bounds<T> evaluate(T x) {
        Bounds<T> sum = T(0, x.precision());

        if (x < 0 || x > 1)
            return sum;

        for (size_t i = 0; i <= _degree; ++i) {
            auto bp = bernsteinBasisPolynomialFor(i, x);
            sum += _coefficients[i] * bp;
        }

        return sum;
    };

    Bounds<T> evaluate(Bounds<T> x) {
        return evaluate(Bounds<T>(x));
    }

  private:
    std::vector<T> computeCoefficients() {
        std::vector<T> coefficients = {};

        for (size_t i = 0; i <= _degree; ++i)
            coefficients.emplace_back(_function(T(Approximation<T>(i / static_cast<float>(_degree), MultiplePrecision(_precision)))));

        return coefficients;
    };

    Bounds<T> bernsteinBasisPolynomialFor(int v, T x) {
        return binomialCoefficients(_degree, v) * pow(x, v) * pow(1 - x, _degree - v);
    }

    const std::function<T(T)>& _function;
    const int _degree;
    std::vector<T> _coefficients;
    const size_t _precision;
};

} // namespace Ariadne

#endif // ARIADNE_BERNSTEIN_POLYNOMIAL_HPP
