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
    BernsteinPolynomial(const std::function<Bounds<T>(Bounds<T>)> &function, DegreeType degree, PR precision, int secantIters = 5) {
        generateCoefficients(function, degree, precision);
        findCriticalPoints(degree, precision, secantIters);
    }

    Bounds<T> evaluate(const Bounds<T> &x) const {
        auto y1 = evaluate_impl(x.lower_raw());
        auto y2 = evaluate_impl(x.upper_raw());

        auto res = Bounds<T>(
            min(y1.lower(), y2.lower()),
            max(y1.upper(), y2.upper()));

        for (const CriticalPoint &criticalPoint: _criticalPoints) {
            if (!models(x, criticalPoint.xPosition))
                continue;

            res = Bounds<T>(
                min(res.lower(), criticalPoint.value.lower()),
                max(res.upper(), criticalPoint.value.upper()));
        }

        return res;
    }

    Bounds<T> evaluateDerivative(const Bounds<T> &x) const {
        auto y1 = evaluate_deriv_impl(x.lower_raw());
        auto y2 = evaluate_deriv_impl(x.upper_raw());

        auto res = Bounds<T>(
            min(y1.lower(), y2.lower()),
            max(y1.upper(), y2.upper()));
        return res;
    }

    Bounds<T> DeCasteljau(const Bounds<T> &x) const {
        std::vector<Bounds<T>> beta = std::vector<Bounds<T>>(_coefficients);

        int n = beta.size();
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < (n - i) - 1; j++) {
                beta[j] = beta[j] * (1 - x) + beta[j + 1] * x;
            }
        }
        return beta[0];
    }

    Bounds<T> operator()(const Bounds<T> &x) const { return evaluate(x); }

  protected:
    Bounds<T> evaluate_impl(const T &x) const {
        if ((x == 1).repr() >= LogicalValue::LIKELY)
            return *(_coefficients.end() - 1);
        else if ((x == 0).repr() >= LogicalValue::LIKELY)
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

    Bounds<T> evaluate_deriv_impl(const T &x) const {
        if ((x == 1).repr() >= LogicalValue::LIKELY)
            return _derivativeEndpoints[1];
        else if ((x == 0).repr() >= LogicalValue::LIKELY)
            return _derivativeEndpoints[0];

        // Derivative of the bernstein is just the sum of basis derivatives (Sum rule)
        // b_(k,n)'(x) = (nCk) * x^(k-1) * (1-x)^(n-k)*  (k-nx)
        auto zero = Bounds<T>(x.precision());
        auto sum = zero;
        Nat degree = _coefficients.size() - 1;

        auto OneMinX = 1 - x;

        auto xPow = 1 / x;
        auto xMinPow = pow(OneMinX, degree - 1);

        auto nx = (degree * x);

        for (size_t i = 0; i <= degree; ++i) {
            auto kMinusNX = i - nx;
            auto bp = xPow * xMinPow * kMinusNX;
            sum = fma(bp, _coefficients[i], sum);

            xPow *= x;
            xMinPow /= OneMinX;
        }

        return sum;
    }

    void generateCoefficients(const std::function<Bounds<T>(Bounds<T>)> &function, DegreeType degree, PR precision) {
        _coefficients.clear();
        _coefficients.reserve(degree + 1);

        Bounds<T> denominator = T(1, precision) / T(degree, precision);
        for (size_t i = 0; i <= degree; ++i) {
            int v2 = (i * 2 >= degree) ? (degree - i) : i;
            auto x = i * denominator;
            auto res = function(x) * binomialCoefficients(degree, v2);

            _coefficients.emplace_back(res);
        }

        {
            auto b = binomialCoefficients(degree, 1);

            _derivativeEndpoints = {
                (function(denominator) - function(0 * denominator)) * b,
                (function((degree - 1) * denominator) - function(degree * denominator)) * b,
            };
        }
    }

    void findCriticalPoints(DegreeType degree, PR precision, int secantIterations = 5) {
        _criticalPoints.clear();

        Bounds<T> denominator = T(1, precision) / T(degree, precision);

        Bounds<T> x = Bounds<T>(LowerBound<T>(0, precision), UpperBound(denominator));

        for (size_t i = 1; i <= degree; ++i, x += denominator) {
            Bounds<T> criticalPoint = secantMethod(x, secantIterations);

            T critX = criticalPoint.value_raw();

            if (!models(x, critX))
                continue;

            Bounds<T> value = this->evaluate_impl(critX);

            _criticalPoints.emplace_back(CriticalPoint(critX, value));
        }
    }

    void findCriticalPoints(DegreeType degree, const PositiveUpperBound<T> &targetEpsilon) {
        _criticalPoints.clear();

        Bounds<T> denominator = T(1, targetEpsilon.precision()) / T(degree, targetEpsilon.precision());

        Bounds<T> x = Bounds<T>(LowerBound<T>(0, targetEpsilon.precision()), UpperBound(denominator));

        for (size_t i = 1; i <= degree; ++i) {
            Bounds<T> criticalPoint = secantMethod(x, targetEpsilon);

            T critX = criticalPoint.value_raw();

            if (!models(x, critX))
                continue;

            Bounds<T> value = this->evaluate_impl(critX);

            _criticalPoints.emplace_back(CriticalPoint(critX, value));

            x += denominator;
        }
    }

    static Bounds<T> bernsteinBasisPolynomialFor(int v, int n, const Bounds<T> &x) {
        return pow(x.value(), v) * pow(1 - x.value(), n - v);
    }

    Bounds<T> secantMethod(const Bounds<T> &x, const PositiveUpperBound<T> &targetEpsilon) const {
        T s[]{
            x.lower_raw(),
            x.upper_raw()};

        while ((mag(s[1] - s[0]) > targetEpsilon).repr() >= LogicalValue::LIKELY)
            if (!secantMethod_impl(s))
                break;
        ;

        auto mini = min(s[0], s[1]);
        auto maxi = max(s[0], s[1]);

        return Bounds<T>(LowerBound<T>(mini), UpperBound<T>(maxi));
    }

    Bounds<T> secantMethod(const Bounds<T> &x, int iterations = 1) const {
        T s[]{
            x.lower_raw(),
            x.upper_raw()};

        for (int i = 0; i < iterations; ++i)
            if (!secantMethod_impl(s))
                break;

        auto mini = min(s[0], s[1]);
        auto maxi = max(s[0], s[1]);

        return Bounds<T>(LowerBound<T>(mini), UpperBound<T>(maxi));
    }

    bool secantMethod_impl(T (&x)[2]) const {
        auto left = x[0];
        auto right = x[1];

        auto rightDeriv = this->evaluate_deriv_impl(right);
        auto leftDeriv = this->evaluate_deriv_impl(left);

        auto res = right - rightDeriv * ((right - left) / (rightDeriv - leftDeriv));

        if (is_nan(res.value()))
            return false;

        x[0] = x[1];
        x[1] = res.value();
        return true;
    }

    std::vector<Bounds<T>> _coefficients{};
    std::vector<Bounds<T>> _derivativeEndpoints;

    struct CriticalPoint {
        CriticalPoint(T xPos, Bounds<T> val) : xPosition(xPos), value(val) {}

        const T xPosition;
        const Bounds<T> value;
    };

    std::vector<CriticalPoint> _criticalPoints{};
};
} // namespace Ariadne

#endif // ARIADNE_BERNSTEIN_POLYNOMIAL_HPP
