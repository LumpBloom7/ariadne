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
class BernsteinPolynomial_impl : virtual public IBernsteinPolynomial<T>,
                                 virtual protected IBernsteinPolynomialBase {
    using RND = T::RoundingModeType;
    using PR = T::PrecisionType;

  public:
    const PR _precision;
    const DegreeType _degree;

    BernsteinPolynomial_impl(std::vector<Bounds<T>> coefficients) : _coefficients{coefficients} {}
    BernsteinPolynomial_impl(const std::function<Bounds<T>(Bounds<T>)> &function, DegreeType degree, PR precision, int secantIters = 5)
        : _precision(precision), _degree(degree), degreeReciprocal(rec(T(degree, precision))), zero(Bounds<T>(precision)) {
        generateCoefficients(function);
        findCriticalPoints(secantIters);
    }

    virtual Bounds<T> evaluate(const Bounds<T> &x) const override {
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

    virtual std::shared_ptr<IBernsteinPolynomial<T>> asSharedPtr() const override {
        return std::make_shared<BernsteinPolynomial_impl<T>>(*this);
    }

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
        auto xMinPow = pow(OneMinX, _degree);

        for (size_t i = 0; i <= _degree; ++i) {
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

        {
            _derivativeEndpoints = {
                (function(degreeReciprocal) - function(zero)) * _degree,
                (function((_degree - 1) * degreeReciprocal) - function(_degree * degreeReciprocal)) * _degree,
            };
        }
    }

    void findCriticalPoints(int secantIterations = 5) {
        _criticalPoints.clear();

        Bounds<T> x = Bounds<T>(zero.lower_raw(), degreeReciprocal.upper_raw());

        for (size_t i = 1; i <= _degree; ++i, x += degreeReciprocal) {
            Bounds<T> criticalPoint = secantMethod(x, secantIterations);

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

        for (size_t i = 1; i <= _degree; ++i) {
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

template<typename T>
class BernsteinPolynomial : public IBernsteinPolynomialPtr<T> {
  public:
    using RND = T::RoundingModeType;
    using PR = T::PrecisionType;

    template<typename... TArgs>
    BernsteinPolynomial(TArgs... args) {
        (this->_ptr) = _ptr2 = std::make_shared<BernsteinPolynomial_impl<T>>(args...);
    }

  private:
    std::shared_ptr<BernsteinPolynomial_impl<T>> _ptr2{nullptr};
};

} // namespace Ariadne

#endif // ARIADNE_BERNSTEIN_POLYNOMIAL_HPP
