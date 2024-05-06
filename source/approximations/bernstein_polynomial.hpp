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
        auto y1 = evaluate_impl(x.lower_raw());
        auto y2 = evaluate_impl(x.upper_raw());

        auto res = Bounds<T>(
            min(y1.lower(), y2.lower()),
            max(y1.upper(), y2.upper()));
        return res;
    }

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
    Bounds<T> evaluate_impl(const T& x) const {
        if ((x >= 1).repr() >= LogicalValue::LIKELY)
            return *(_coefficients.end() - 1);
        else if ((x <= 0).repr() >= LogicalValue::LIKELY)
            return _coefficients[0];

        auto zero = Bounds<T>(x.precision());
        auto sum = zero;

        Nat degree = _coefficients.size() - 1;

        auto OneMinX = 1 - x;

        auto xPow = pow(x, 0);
        auto xMinPow = pow(OneMinX, degree);

        for (size_t i = 0; i <= degree; ++i) {
            auto bp = xPow * xMinPow;
            sum = fma(bp, _coefficients[i], sum);

            xPow *= x;
            xMinPow /= OneMinX;
        }

        return sum;
    }

    void generateCoefficients(const std::function<Bounds<T>(Bounds<T>)>& function, DegreeType degree, PR precision) {
        _coefficients.clear();
        _coefficients.reserve(degree + 1);

        Bounds<T> denominator = T(1, precision) / T(degree, precision);
        for (size_t i = 0; i <= degree; ++i) {
            int v2 = (i * 2 >= degree) ? (degree - i) : i;
            auto x = i * denominator;
            auto res = function(x) * binomialCoefficients(degree, v2);

            _coefficients.emplace_back(res);
        }
    }

    static Bounds<T> bernsteinBasisPolynomialFor(int v, int n, const Bounds<T>& x) {
        return pow(x.value(), v) * pow(1 - x.value(), n - v);
    }

    std::vector<Bounds<T>> _coefficients{};

    static bool straddles(const Bounds<T>& range, const FloatMP& point) {
        return (range.upper_raw() >= point && range.lower_raw() <= point).repr() >= LogicalValue::LIKELY;
    }
};
} // namespace Ariadne

#endif // ARIADNE_BERNSTEIN_POLYNOMIAL_HPP
