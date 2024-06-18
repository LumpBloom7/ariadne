#ifndef ARIADNE_BERNSTEIN_POLYNOMIAL_HPP
#define ARIADNE_BERNSTEIN_POLYNOMIAL_HPP

#include <cmath>
#include <functional>
#include <stack>
#include <vector>

#include "approximations/polynomial_approximation_interface.hpp"
#include "numeric/builtin.hpp"
#include "numeric/float_bounds.hpp"
#include "numeric/integer.hpp"
#include "numeric/real.hpp"
#include "utility/cache.hpp"
#include "utility/factorials.hpp"
#include "utility/hash.hpp"
#include "utility/hash_gmp.hpp"
#include "utility/hash_tuple.hpp"

namespace Ariadne {

class IBernsteinPolynomialBase {
  protected:
    static SimpleCache<Integer, Nat64, Nat64> binomialCoefficients;
};

template<typename T>
class BernsteinPolynomial : virtual public IPolynomialApproximation<T>,
                            virtual protected IBernsteinPolynomialBase {
    using RND = T::RoundingModeType;
    using PR = T::PrecisionType;

  public:
    PR _precision;
    DegreeType _degree;

    BernsteinPolynomial(std::vector<Bounds<T>> coefficients)
        : _coefficients{coefficients}, _precision(coefficients[0].precision()), _degree(coefficients.size()),
          degreeReciprocal(rec(T(_degree, _precision))), zero(Bounds<T>(_precision)) {}
    BernsteinPolynomial(const std::function<Bounds<T>(Bounds<T>)> &function, DegreeType degree, PR precision)
        : _precision(precision), _degree(degree), degreeReciprocal(rec(T(degree, precision))),
          zero(Bounds<T>(precision)) {
        generateCoefficients(function);
        computeDerivEndpoints();
    }

#ifndef _OPENMP
    Bounds<T> range(const Interval<T> &x, Nat subIntervals = 1) const final {
        auto stepSize = (x.upper_raw() - x.lower_raw()) / subIntervals;

        auto subinterval = Bounds<T>(x.lower_raw(), (x.lower_raw() + stepSize).upper_raw());

        auto y1 = evaluate_impl(fromInterval(subinterval));

        auto mini = y1.lower_raw();
        auto maxi = y1.upper_raw();

        for (int i = 1; i < subIntervals; ++i) {
            subinterval += stepSize;
            auto res = evaluate_impl(subinterval);

            mini = min(mini, res.lower_raw()), maxi = max(maxi, res.upper_raw());
        }

        return Bounds<T>(mini, maxi);
    }
#else
    Bounds<T> range(const Interval<T> &x, Nat subIntervals = 1) const final {
        // Skip OpenMP overheads
        if (subIntervals == 1)
            return evaluate_impl(this->fromInterval(x));

        auto stepSize = (x._u - x._l) / subIntervals;
        auto interval = Bounds<T>(x._l, (x._l + stepSize).upper_raw());

        std::vector<Bounds<T>> results{};

        for (int i = 0; i < subIntervals; ++i)
            results.emplace_back(zero);

#pragma omp parallel for
        for (int i = 0; i < subIntervals; ++i) {
            auto subinterval = interval + (stepSize * i);

            results[i] = evaluate_impl(interval);
        }

        auto mini = results[0].lower_raw();
        auto maxi = results[0].upper_raw();
        for (int i = 1; i < subIntervals; ++i) {
            auto &res = results[i];
            mini = min(mini, res.lower_raw());
            maxi = max(maxi, res.upper_raw());
        }

        return Bounds<T>(mini, maxi);
    }
#endif

    Bounds<T> evaluate(const Bounds<T> &x) const final { return evaluate_impl(x); }

    Bounds<T> evaluateDerivative(const Bounds<T> &x) const final {
        auto y1 = evaluate_deriv_impl(x.lower_raw());
        auto y2 = evaluate_deriv_impl(x.upper_raw());

        auto res = Bounds<T>(min(y1.lower_raw(), y2.lower_raw()), max(y1.upper(), y2.upper()));

        return res;
    }

    DegreeType degree() const final { return _degree; }

    PR precision() const final { return _precision; }

    std::vector<T> coefficients() const final { return _coefficients; }

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
            _coefficients[i] = cast_exact(Approximation<T>(_coefficients[i] + other._coefficients[i]));

        return BernsteinPolynomial<T>(coefficients);
    }

    BernsteinPolynomial<T> &operator+=(const BernsteinPolynomial<T> &other) {
        if (other._degree != _degree)
            throw std::runtime_error("Degree of both polynomials do not match");

        for (int i = 0; i < _degree; ++i)
            _coefficients[i] = cast_exact(Approximation<T>(_coefficients[i] + other._coefficients[i]));

        computeDerivEndpoints();
        return *this;
    }

    BernsteinPolynomial<T> &operator-=(const BernsteinPolynomial<T> &other) {
        if (other._degree != _degree)
            throw std::runtime_error("Degree of both polynomials do not match");

        for (int i = 0; i < _degree; ++i)
            _coefficients[i] -= other._coefficients[i];

        computeDerivEndpoints();
        return *this;
    }

  protected:
    Bounds<T> evaluate_impl(const Bounds<T> &x) const {
        auto sum = zero;

        auto OneMinX = 1 - x;

        auto xPow = pow(x, 0);
        auto OneMinXPow = pow(OneMinX, 0);

        std::stack<Bounds<T>> oneMinXPow = {};

        // Computing (1-x)^n, then dividing by (1-x) to reduce the power is not
        // ideal due to potential (1/0) error Instead we do it backwards, and store
        // the result in a stack (Multiplying by 0 is well defined) This incurs a
        // slight memory overhead, but the benefit of removing pow(1-x, i) is worth
        // it.
        for (int i = 0; i <= _degree; ++i) {
            oneMinXPow.emplace(OneMinXPow);
            OneMinXPow *= OneMinX;
        }

        for (size_t i = 0; i <= _degree; ++i) {
            // xPow is multiplied incrementally, this is fine a 0^x is always valid
            // I want to do the same for the 1-x part, but 0^-1 is not valid, and the
            // default pow function handles is just fine
            auto bp = xPow * oneMinXPow.top();

            oneMinXPow.pop();

            sum = fma(bp, _coefficients[i], sum);

            xPow *= x;
        }
        return sum;
    }

    Bounds<T> evaluate_deriv_impl(const T &x) const {
        if ((x == 1).repr() >= LogicalValue::LIKELY)
            return _derivativeEndpoints[1];
        else if ((x == 0).repr() >= LogicalValue::LIKELY)
            return _derivativeEndpoints[0];

        // Derivative of the bernstein is just the sum of basis derivatives (Sum
        // rule) b_(k,n)'(x) = (nCk) * x^(k-1) * (1-x)^(n-k)*  (k-nx)
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

            auto approximation = Approximation<T>(res);
            auto exact = cast_exact(approximation);

            _coefficients.emplace_back(exact);
            x += degreeReciprocal;
        }
    }

    void computeDerivEndpoints() {
        _derivativeEndpoints = {(_coefficients[1] - _coefficients[0]) * _degree,
                                (_coefficients[_degree - 1] - _coefficients[_degree]) * _degree};
    }

    static Bounds<T> bernsteinBasisPolynomialFor(int v, int n, const Bounds<T> &x) {
        return pow(x, v) * pow(1 - x, n - v);
    }

    std::vector<T> _coefficients{};
    std::vector<Bounds<T>> _derivativeEndpoints;
    Bounds<T> degreeReciprocal;
    Bounds<T> zero;
};

} // namespace Ariadne

#endif // ARIADNE_BERNSTEIN_POLYNOMIAL_HPP
