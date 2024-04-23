#ifndef ARIADNE_BERNSTEIN_POLYNOMIAL_HPP
#define ARIADNE_BERNSTEIN_POLYNOMIAL_HPP

#include <cmath>
#include <functional>
#include <iostream>
#include <optional>

#include "numeric/float_approximation.hpp"
#include "numeric/float_bounds.hpp"
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
    // BernsteinPolynomial(std::vector<Bounds<T>> coefficients){}

    BernsteinPolynomial(std::function<Bounds<T>(Bounds<T>)> function, int degree) : _function{function}, _degree{degree} {}

    PositiveUpperBound<T> maximumError() {
        auto maximum = _errorBounds[0];

        for (int i = 1; i < _errorBounds.size(); ++i) {
            auto bounds = _errorBounds[i];

            maximum = max(bounds, maximum);
        }

        return maximum;
    }
    PositiveUpperBound<T> maximumErrorAt(Bounds<T> x) {
        T denominator = T(Approximation<T>(1 / static_cast<float>(_degree), x.precision()));

        auto maximum = PositiveUpperBound<T>(0.0_x, x.precision());
 
        for (size_t i = 0; i < _degree; ++i) {
            auto xl = i * denominator;
            auto xr = (i + 1) * denominator;

            if ((xl > x.upper() || xr < x.lower()).repr() >= LogicalValue::LIKELY)
                continue;

            maximum = max(maximum, _errorBounds[i]);
        }

        return maximum;
    }

    Bounds<T> evaluate(Bounds<T> x) {
        // This is out of evaluation domain, return 0.
        // We could remove this with little consequence to correctness, as evaluating a BernsteinPolynomial like this would be undefined anyways.
        // Doing this allows us to skip computations tExactDoublehough.
        if (Detail::possibly((x < 0 || x > 1).repr()))
            return T(0, x.precision());

        validateCoefficients(x.precision());

        auto sum = Bounds<T>(x.precision());

        for (size_t i = 0; i <= _degree; ++i) {
            auto bp = bernsteinBasisPolynomialFor(i, x);
            sum += _coefficients[i] * bp;
        }

        return sum;
    };

  private:
    void validateCoefficients(PR precisionType) {
        // Coefficients are sufficient for the desired evaluation precision. Avoid recomputation.
        if (!_coefficients.empty() && _coefficients[0].precision() >= precisionType)
            return;

        _coefficients = {};
        _errorBounds = {};

        T denominator = T(Approximation<T>(1 / static_cast<float>(_degree), precisionType));
        _coefficients.emplace_back(_function(Bounds<T>(0, precisionType)));

        for (size_t i = 1; i <= _degree; ++i) {
            auto x = i * denominator;
            auto res = _function(i * denominator);

            auto interval = Bounds<T>(LowerBound<T>((i - 1) / denominator), UpperBound<T>(x));

            _errorBounds.emplace_back(mag(_function(interval)));
            _coefficients.emplace_back(res);
        }
    };

    Bounds<T> bernsteinBasisPolynomialFor(int v, Bounds<T> x) {
        return binomialCoefficients(_degree, v) * pow(x, v) * pow(1 - x, _degree - v);
    }

    std::function<Bounds<T>(Bounds<T>)> _function;
    const int _degree;
    std::vector<Bounds<T>> _coefficients;
    std::vector<PositiveUpperBound<T>> _errorBounds;
};

// TODO: We can't approximate easily using Reals, as it would require several computations of the original functions, which is not ideal
// I want to cache the coefficient as Either a validated real, or a dyadic bounds. But Reals don't support math operations with those types straightforward.

/*
template<>
class BernsteinPolynomial<Real> : protected BernsteinPolynomialBase {
  public:
    BernsteinPolynomial(const std::function<Real(Real)>& function, int degree, size_t max_precision) : _function{function}, _degree{degree}, _precision{max_precision} {
        _coefficients = computeCoefficients();
    }

    DyadicBounds evaluate(Real x) {
        DyadicBounds sum = DyadicBounds(0);

        for (size_t i = 0; i <= _degree; ++i) {
            auto bp = bernsteinBasisPolynomialFor(i, x);
            sum += _coefficients[i] * bp;
        }

        return sum;
    };

  private:
    std::vector<Bounds<Dyadic>> computeCoefficients() {
        std::vector<Bounds<Dyadic>> coefficients = {};

        auto D = Real(_degree);

        for (size_t i = 0; i <= _degree; ++i) {
            coefficients.emplace_back(_function(Real(i) / D).compute(Effort(_precision)));
        }

        return coefficients;
    };

    Real bernsteinBasisPolynomialFor(int v, Real x) {
        return binomialCoefficients(_degree, v) * pow(x, v) * pow(1 - x, _degree - v);
    }

    const std::function<Real(Real)>& _function;
    const int _degree;
    std::vector<Bounds<Dyadic>> _coefficients;
    const size_t _precision;
}; */
} // namespace Ariadne

#endif // ARIADNE_BERNSTEIN_POLYNOMIAL_HPP
