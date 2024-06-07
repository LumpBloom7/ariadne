#ifndef ARIADNE_BERNSTEIN_POLYNOMIAL_HPP
#define ARIADNE_BERNSTEIN_POLYNOMIAL_HPP

#include <cmath>
#include <functional>
#include <iostream>
#include <optional>

#include "approximations/bernstein_polynomial_shared.hpp"
#include "numeric/float_approximation.hpp"
#include "numeric/float_bounds.hpp"
#include "numeric/float_error.hpp"
#include "numeric/floatmp.hpp"
#include "numeric/integer.hpp"
#include "numeric/real.hpp"
#include "numeric/validated_real.hpp"
#include "utility/hash.hpp"
#include "utility/hash_numeric.hpp"
#include "utility/hash_tuple.hpp"
#include "utility/standard.hpp"
namespace Ariadne {

template<typename T>
class BernsteinPolynomial : virtual public IPolynomialApproximation<T>,
                            virtual protected IBernsteinPolynomialBase {
    using RND = T::RoundingModeType;
    using PR = T::PrecisionType;

  public:
    PR _precision;
    DegreeType _degree;

    BernsteinPolynomial(std::vector<Bounds<T>> coefficients) : _coefficients{coefficients},
                                                               _precision(coefficients[0].precision()), _degree(coefficients.size()), degreeReciprocal(rec(T(_degree, _precision))), zero(Bounds<T>(_precision)) {}
    BernsteinPolynomial(const std::function<Bounds<T>(Bounds<T>)> &function, DegreeType degree, PR precision)
        : _precision(precision), _degree(degree), degreeReciprocal(rec(T(degree, precision))), zero(Bounds<T>(precision)) {
        generateCoefficients(function);
        computeDerivEndpoints();
    }

    /*     virtual Bounds<T> evaluate(const Bounds<T> &x) const override {
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
     */
    virtual Bounds<T> evaluate(const Bounds<T> &x, int subIntervals = 1) const override {
        auto stepSize = (x.upper_raw() - x.lower_raw()) / subIntervals;

        auto subinterval = Bounds<T>(x.lower_raw(), (x.lower_raw() + stepSize).upper_raw());

        auto y1 = evaluate_impl(subinterval);

        auto mini = y1.lower_raw();
        auto maxi = y1.upper_raw();

        for (int i = 1; i < subIntervals; ++i) {
            subinterval += stepSize;
            auto res = evaluate_impl(subinterval);

            mini = min(mini, res.lower_raw()),
            maxi = max(maxi, res.upper_raw());
        }

        return Bounds<T>(mini, maxi);
    }

    virtual Bounds<T> evaluateRaw(const Bounds<T> &x) const override {
        return evaluate_impl(x);
    }

    virtual Bounds<T> evaluateDerivative(const Bounds<T> &x) const override {
        auto y1 = evaluate_deriv_impl(x.lower_raw());
        auto y2 = evaluate_deriv_impl(x.upper_raw());

        auto res = Bounds<T>(
            min(y1.lower_raw(), y2.lower_raw()),
            max(y1.upper(), y2.upper()));

        return res;
    }

    virtual DegreeType degree() const override {
        return _degree;
    };

    virtual PR precision() const override {
        return _precision;
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

    BernsteinPolynomial<T> operator+(const BernsteinPolynomial<T> &other) const {
        if (other._degree != _degree)
            throw std::runtime_error("Degree of both polynomials do not match");

        auto coefficients = _coefficients;

        for (int i = 0; i < _degree; ++i)
            coefficients[i] += other._coefficients[i];

        return BernsteinPolynomial<T>(coefficients);
    }

    BernsteinPolynomial<T> operator-(const BernsteinPolynomial<T> &other) const {
        if (other._degree != _degree)
            throw std::runtime_error("Degree of both polynomials do not match");

        auto coefficients = _coefficients;

        for (int i = 0; i < _degree; ++i)
            coefficients[i] -= other._coefficients[i];

        return BernsteinPolynomial<T>(coefficients);
    }

    BernsteinPolynomial<T> &operator+=(const BernsteinPolynomial<T> &other) {
        if (other._degree != _degree)
            throw std::runtime_error("Degree of both polynomials do not match");

        for (int i = 0; i < _degree; ++i)
            _coefficients[i] += other._coefficients[i];

        computeDerivEndpoints();
        findCriticalPoints();
        return *this;
    }

    BernsteinPolynomial<T> &operator-=(const BernsteinPolynomial<T> &other) {
        if (other._degree != _degree)
            throw std::runtime_error("Degree of both polynomials do not match");

        for (int i = 0; i < _degree; ++i)
            _coefficients[i] -= other._coefficients[i];

        computeDerivEndpoints();
        findCriticalPoints();
        return *this;
    }

  protected:
    Bounds<T> evaluate_impl(const Bounds<T> &x) const {
        auto sum = zero;
        for (size_t i = 0; i <= _degree; ++i) {
            auto bp = bernsteinBasisPolynomialFor(i, _degree, x);
            sum = fma(bp, _coefficients[i], sum);
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
        auto xMinPow = pow(OneMinX, _degree - 1);

        auto nx = (_degree * x);

        for (size_t i = 0; i <= _degree; ++i) {
            auto kMinusNX = i - nx;
            auto bp = xPow * xMinPow * kMinusNX;
            sum = fma(bp, _coefficients[i], sum);

            xPow *= x;
            xMinPow *= oneMinXRec;
        }

        return sum;
    }

    void generateCoefficients(const std::function<Bounds<T>(Bounds<T>)> &function) {
        _coefficients.clear();
        _coefficients.reserve(_degree + 1);

        auto x = zero;
        for (size_t i = 0; i <= _degree; ++i) {
            int v2 = (i * 2 >= _degree) ? (_degree - i) : i;
            auto res = function(x) * binomialCoefficients(_degree, v2);

            _coefficients.emplace_back(res);
            x += degreeReciprocal;
        }
    }

    void computeDerivEndpoints() {
        _derivativeEndpoints = {
            (_coefficients[1] - _coefficients[0]) * _degree,
            (_coefficients[_degree - 1] - _coefficients[_degree]) * _degree};
    }

    void findCriticalPoints(int secantIterations = 5) {
        _criticalPoints.clear();

        if (secantIterations == 0)
            return;

        Bounds<T> x = Bounds<T>(zero.lower_raw(), degreeReciprocal.upper_raw());

        T lastCritX = x.lower_raw();

        for (size_t i = 1; i <= _degree; ++i, x += degreeReciprocal) {
            if (decide(x.upper_raw() < lastCritX)) {
                continue;
            }

            Bounds<T> criticalPoint = secantMethod(x, secantIterations);
            lastCritX = criticalPoint.lower_raw();

            T critX = criticalPoint.value_raw();

            if (!models(x, critX))
                continue;

            Bounds<T> value = this->evaluate_impl(critX);

            _criticalPoints.emplace_back(CriticalPoint(critX, value));
        }
    }

    void findCriticalPoints(const PositiveUpperBound<T> &targetEpsilon) {
        _criticalPoints.clear();

        Bounds<T> x = Bounds<T>(zero.lower_raw(), degreeReciprocal);

        T lastCritX = x.value_raw();

        for (size_t i = 1; i <= _degree; ++i) {
            if (decide(x.upper_raw() < lastCritX))
                continue;

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
        return pow(x, v) * pow(1 - x, n - v);
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

        for (int i = 0; i < iterations; ++i) {
            if (decide(max(s[0], s[1]) <= x.lower_raw()))
                break;

            if (!secantMethod_impl(s, leftval))
                break;
        }

        auto mini = min(s[0], s[1]);
        auto maxi = max(s[0], s[1]);

        return Bounds<T>(mini, maxi);
    }

    bool secantMethod_impl(T (&x)[2], Bounds<T> &leftDeriv) const {
        auto &left = x[0];
        auto &right = x[1];

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
    Bounds<T> zero;

    struct CriticalPoint {
        CriticalPoint(T xPos, Bounds<T> val) : xPosition(xPos), value(val) {}

        const T xPosition;
        const Bounds<T> value;
    };

    std::vector<CriticalPoint> _criticalPoints{};
};

} // namespace Ariadne

#endif // ARIADNE_BERNSTEIN_POLYNOMIAL_HPP
