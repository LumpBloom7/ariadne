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
    const PR precision;
    const DegreeType degree;

    BernsteinPolynomial(std::vector<Bounds<T>> coefficients) : _coefficients{coefficients} {}
    BernsteinPolynomial(const std::function<Bounds<T>(Bounds<T>)> &function, DegreeType degree, PR precision, int secantIters = 5)
        : precision(precision), degree(degree), degreeReciprocal(rec(T(degree, precision))), zero(Bounds<T>(precision)) {
        generateCoefficients(function, degree, precision);
        findCriticalPoints(degree, precision, secantIters);
    }

    Bounds<T> evaluate(const Bounds<T> &x) const {
        auto y1 = evaluate_impl(x.lower_raw());
        auto y2 = evaluate_impl(x.upper_raw());

        auto mini = min(y1.lower_raw(), y2.upper_raw());
        auto maxi = max(y1.upper_raw(), y2.upper_raw());

        for (const CriticalPoint &criticalPoint: _criticalPoints) {
            if (decide(criticalPoint.xPosition > x.upper_raw()))
                break;

            if (criticalPoint.xPosition < x.lower_raw())
                continue;

            mini = min(mini, criticalPoint.value.lower_raw()),
            maxi = max(maxi, criticalPoint.value.upper_raw());
        }

        return Bounds<T>(mini, maxi);
    }

    Bounds<T> evaluateDerivative(const Bounds<T> &x) const {
        auto y1 = evaluate_deriv_impl(x.lower_raw());
        auto y2 = evaluate_deriv_impl(x.upper_raw());

        auto res = Bounds<T>(
            min(y1.lower_raw(), y2.lower_raw()),
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

        auto sum = zero;

        auto OneMinX = 1 - x;
        auto oneMinXRec = rec(OneMinX);

        auto xPow = pow(x, 0);
        auto xMinPow = pow(OneMinX, degree);

        for (size_t i = 0; i <= degree; ++i) {
            auto bp = xPow * xMinPow;
            sum = fma(bp, _coefficients[i], sum);

            xPow *= x;
            xMinPow *= oneMinXRec;
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
        auto sum = zero;

        auto OneMinX = 1 - x;
        auto oneMinXRec = rec(OneMinX);

        auto xPow = 1 / x;
        auto xMinPow = pow(OneMinX, degree - 1);

        auto nx = (degree * x);

        for (size_t i = 0; i <= degree; ++i) {
            auto kMinusNX = i - nx;
            auto bp = xPow * xMinPow * kMinusNX;
            sum = fma(bp, _coefficients[i], sum);

            xPow *= x;
            xMinPow *= oneMinXRec;
        }

        return sum;
    }

    void generateCoefficients(const std::function<Bounds<T>(Bounds<T>)> &function, DegreeType degree, PR precision) {
        _coefficients.clear();
        _coefficients.reserve(degree + 1);

        auto x = zero;
        for (size_t i = 0; i <= degree; ++i) {
            int v2 = (i * 2 >= degree) ? (degree - i) : i;
            auto res = function(x) * binomialCoefficients(degree, v2);

            _coefficients.emplace_back(res);
            x += degreeReciprocal;
        }

        {
            _derivativeEndpoints = {
                (function(degreeReciprocal) - function(zero)) * degree,
                (function((degree - 1) * degreeReciprocal) - function(degree * degreeReciprocal)) * degree,
            };
        }
    }

    void findCriticalPoints(DegreeType degree, PR precision, int secantIterations = 5) {
        _criticalPoints.clear();

        Bounds<T> x = Bounds<T>(zero.lower_raw(), degreeReciprocal.upper_raw());

        for (size_t i = 1; i <= degree; ++i, x += degreeReciprocal) {
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

        Bounds<T> x = Bounds<T>(zero.lower_raw(), degreeReciprocal);

        for (size_t i = 1; i <= degree; ++i) {
            Bounds<T> criticalPoint = secantMethod(x, targetEpsilon);

            T critX = criticalPoint.value_raw();

            if (!models(x, critX))
                continue;

            Bounds<T> value = this->evaluate_impl(critX);

            _criticalPoints.emplace_back(CriticalPoint(critX, value));

            x += degreeReciprocal;
        }
    }

    static Bounds<T> bernsteinBasisPolynomialFor(int v, int n, const Bounds<T> &x) {
        return pow(x.value(), v) * pow(1 - x.value(), n - v);
    }

    const Bounds<T> secantMethod(const Bounds<T> &x, const PositiveUpperBound<T> &targetEpsilon) const {
        T s[]{
            x.lower_raw(),
            x.upper_raw()};

        auto leftval = this->evaluate_deriv_impl(s[0]);

        while ((mag(s[1] - s[0]) > targetEpsilon).repr() >= LogicalValue::LIKELY)
            if (!secantMethod_impl(s, leftval))
                break;
        ;

        auto mini = min(s[0], s[1]);
        auto maxi = max(s[0], s[1]);

        return Bounds<T>(mini, maxi);
    }

    const Bounds<T> secantMethod(const Bounds<T> &x, int iterations = 1) const {
        T s[]{
            x.lower_raw(),
            x.upper_raw()};

        auto leftval = this->evaluate_deriv_impl(s[0]);

        for (int i = 0; i < iterations; ++i)
            if (!secantMethod_impl(s, leftval))
                break;

        auto mini = min(s[0], s[1]);
        auto maxi = max(s[0], s[1]);

        return Bounds<T>(mini, maxi);
    }

    bool secantMethod_impl(T (&x)[2], Bounds<T> &leftDeriv) const {
        auto left = x[0];
        auto right = x[1];

        auto rightDeriv = this->evaluate_deriv_impl(right);

        auto res = (right - rightDeriv * ((right - left) / (rightDeriv - leftDeriv))).value();

        if (is_nan(res))
            return false;

        x[0] = x[1];
        x[1] = res;
        leftDeriv = rightDeriv;
        return true;
    }

    std::vector<Bounds<T>> _coefficients{};
    std::vector<Bounds<T>> _derivativeEndpoints;
    Bounds<T> degreeReciprocal;
    const Bounds<T> zero;

    struct CriticalPoint {
        CriticalPoint(T xPos, Bounds<T> val) : xPosition(xPos), value(val) {}

        const T xPosition;
        const Bounds<T> value;
    };

    std::vector<CriticalPoint> _criticalPoints{};
};
} // namespace Ariadne

#endif // ARIADNE_BERNSTEIN_POLYNOMIAL_HPP
