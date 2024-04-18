#ifndef ARIADNE_BERNSTEIN_POLYNOMIAL_HPP
#define ARIADNE_BERNSTEIN_POLYNOMIAL_HPP

#include <cmath>
#include <functional>
#include <iostream>

#include "numeric/float_approximation.hpp"
#include "numeric/float_bounds.hpp"
#include "numeric/floatmp.hpp"
#include "numeric/integer.hpp"
#include "utility/cache.hpp"
#include "utility/factorials.hpp"
#include "utility/hash.hpp"
#include "utility/hash_numeric.hpp"
#include "utility/hash_tuple.hpp"
#include "utility/standard.hpp"

namespace Ariadne {

class BernsteinPolynomialBase {
  protected:
    static SimpleCache<Integer, Nat64, Nat64> binomialCoefficients;
};

template<typename T>
    requires HasPrecisionType<T> && ConvertibleTo<T, Bounds<T>>
class BernsteinPolynomial : protected BernsteinPolynomialBase {
  public:
    BernsteinPolynomial(const std::function<Bounds<T>(Bounds<T>)>& function, int degree, size_t max_precision) : _function{function}, _degree{degree}, _precision{max_precision} {
        _coefficients = computeCoefficients();
    }

    Bounds<T> evaluate(Bounds<T> x) {
        Bounds<T> sum = T(0, x.precision());

        if (Detail::possibly((x < 0 ||x > 1).repr()))
            return sum;

        for (size_t i = 0; i <= _degree; ++i) {
            auto bp = bernsteinBasisPolynomialFor(i, x);
            sum += _coefficients[i] * bp;
        }

        return sum;
    };

  private:
    std::vector<Bounds<T>> computeCoefficients() {
        std::vector<Bounds<T>> coefficients = {};

        for (size_t i = 0; i <= _degree; ++i)
            coefficients.emplace_back(_function(T(Approximation<T>(i / static_cast<float>(_degree), MultiplePrecision(_precision)))));

        return coefficients;
    };

    Bounds<T> bernsteinBasisPolynomialFor(int v, Bounds<T> x) {
        return binomialCoefficients(_degree, v) * pow(x, v) * pow(1 - x, _degree - v);
    }

    const std::function<Bounds<T>(Bounds<T>)>& _function;
    const int _degree;
    std::vector<Bounds<T>> _coefficients;
    const size_t _precision;
};

} // namespace Ariadne

#endif // ARIADNE_BERNSTEIN_POLYNOMIAL_HPP
