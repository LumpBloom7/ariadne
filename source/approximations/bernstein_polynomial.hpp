#ifndef ARIADNE_BERNSTEIN_POLYNOMIAL_HPP
#define ARIADNE_BERNSTEIN_POLYNOMIAL_HPP

#include <cmath>
#include <functional>
#include <iostream>
#include <optional>

#include "numeric/float_approximation.hpp"
#include "numeric/float_bounds.hpp"
#include "numeric/float_error.hpp"
#include "numeric/floatmp.hpp"
#include "numeric/integer.hpp"
#include "numeric/real.hpp"
#include "numeric/validated_real.hpp"
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
class BernsteinPolynomial : protected BernsteinPolynomialBase {
    typedef typename T::RoundingModeType RND;
    typedef typename T::PrecisionType PR;

  public:
    BernsteinPolynomial(std::vector<Bounds<T>> coefficients) : _coefficients{coefficients} {}
    BernsteinPolynomial(const std::function<Bounds<T>(Bounds<T>)>& function, DegreeType degree, PR precision) {
        generateCoefficients(function, degree, precision);
    }

    Bounds<T> evaluate(const Bounds<T>& x) const {
        auto sum = Bounds<T>(x.precision());

        auto degree = _coefficients.size() - 1;
        for (size_t i = 0; i <= degree; ++i) {
            auto bp = bernsteinBasisPolynomialFor(i, degree, x);
            sum = fma(bp, _coefficients[i], sum);
        }

        return sum;
    };

    Bounds<T> DeCasteljau(const Bounds<T>& x) const {
        std::vector<Bounds<T>> beta = std::vector<Bounds<T>>(_coefficients);

        int n = beta.size();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < (n - i) - 1; j++) {
                beta[j] = beta[j] * (1 - x) + beta[j + 1] * x;
            }
        }
        return beta[0];
    }

    Bounds<T> operator()(const Bounds<T>& x) const { return evaluate(x); }

  protected:
    void generateCoefficients(const std::function<Bounds<T>(Bounds<T>)>& function, DegreeType degree, PR precision) {
        _coefficients.clear();
        _coefficients.reserve(degree + 1);

        Bounds<T> denominator = T(1, precision) / T(degree, precision);
        for (size_t i = 0; i <= degree; ++i) {
            auto x = i * denominator;
            auto res = function(x);

            _coefficients.emplace_back(res);
        }
    }

    static Bounds<T> bernsteinBasisPolynomialFor(int v, int n, const Bounds<T>& x) {
        int v2 = (v * 2 >= n) ? (n - v) : v; // Used to avoid computing extra redundant values in binomial cache

        return binomialCoefficients(n, v2) * pow(x, v) * pow(1 - x, n - v);
    }
    std::vector<Bounds<T>> _coefficients{};
};
} // namespace Ariadne

#endif // ARIADNE_BERNSTEIN_POLYNOMIAL_HPP
