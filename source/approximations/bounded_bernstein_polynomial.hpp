#ifndef ARIADNE_BOUNDED_BERNSTEIN_POLYNOMIAL_HPP
#define ARIADNE_BOUNDED_BERNSTEIN_POLYNOMIAL_HPP

#include <cmath>
#include <functional>
#include <iostream>
#include <optional>

#include "approximations/bernstein_polynomial.hpp"

namespace Ariadne {
template<typename T>
class BoundedBernsteinPolynomial : public BernsteinPolynomial<T> {
    typedef typename T::RoundingModeType RND;
    typedef typename T::PrecisionType PR;

  public:
    BoundedBernsteinPolynomial(const std::function<Bounds<T>(Bounds<T>)>& function, DegreeType degree, PR precision)
        : BernsteinPolynomial<T>(function, degree, precision) {
        computeErrorBounds(function, PositiveUpperBound<T>(T::inf(precision)));
    }

    BoundedBernsteinPolynomial(const std::function<Bounds<T>(Bounds<T>)>& function, const PositiveUpperBound<T>& targetEpsilon)
        : BernsteinPolynomial<T>({}) {
        DegreeType degree = 1;

        while (!earlyBoundsTest(function, degree, targetEpsilon)) {
            degree *= 2;
        }

        std::cout << degree << std::endl;

        do {
            this->generateCoefficients(function, degree, targetEpsilon.precision());
            degree *= 2;
        } while (!firstPassTest(function, targetEpsilon) || !computeErrorBounds(function, targetEpsilon));
    }

    PositiveUpperBound<T> maximumError() const {
        auto maximum = _errorBounds[0];

        for (int i = 1; i < _errorBounds.size(); ++i) {
            maximum = max(maximum, _errorBounds[i]);
        }

        return maximum;
    }

    PositiveUpperBound<T> maximumErrorAt(const Bounds<T>& x) const {
        auto degree = (this->_coefficients).size() - 1;
        auto denominator = invert(T(degree, x.precision()));
        auto maximum = PositiveUpperBound<T>(x.precision());

        for (size_t i = 0; i <= degree; ++i) {
            auto xl = i * denominator;
            auto xr = (i + 1) * denominator;

            if ((xl > x.upper() || xr < x.lower()).repr() >= LogicalValue::LIKELY)
                continue;

            maximum = max(maximum, _errorBounds[i]);
        }
        std::cout << std::endl;

        return maximum;
    }

  private:
    bool earlyBoundsTest(const std::function<Bounds<T>(Bounds<T>)>& function, DegreeType degree, const PositiveUpperBound<T>& targetEpsilon) {
        auto denominator = invert(T(degree, targetEpsilon.precision()));
        for (int i = 1; i <= degree; ++i) {
            auto interval = Bounds<T>(LowerBound<T>((i - 1) * denominator), UpperBound<T>(i * denominator));

            auto range = function(interval);
            auto error = mag(range.upper().raw() - range.lower().raw());

            std::cout << interval << std::endl;
            std::cout << range << std::endl;
            std::cout << error << std::endl
                      << degree << std::endl
                      << std::endl;

            if ((error > targetEpsilon).repr() >= LogicalValue::INDETERMINATE)
                return false;
        }

        return true;
    }

    bool firstPassTest(const std::function<Bounds<T>(Bounds<T>)>& function, const PositiveUpperBound<T>& targetEpsilon) {
        auto degree = (this->_coefficients).size() - 1;
        auto denominator = invert(T(degree, targetEpsilon.precision()));

        for (int i = 1; i < degree; ++i) {
            auto x = i * denominator;

            auto actual = function(x);
            auto predicted = this->evaluate(x);

            auto maxError = mag(actual - predicted);

/*             std::cout << x << std::endl;
            std::cout << actual << std::endl;
            std::cout << predicted << std::endl;
            std::cout << maxError << std::endl
                      << degree << std::endl
                      << std::endl; */

            if ((maxError > targetEpsilon).repr() >= LogicalValue::INDETERMINATE) {
                return false;
            }
        }

        return true;
    }

    bool computeErrorBounds(const std::function<Bounds<T>(Bounds<T>)>& function, const PositiveUpperBound<T>& targetEpsilon) {
        auto degree = (this->_coefficients).size() - 1;
        _errorBounds.clear();
        _errorBounds.reserve(degree);
        auto denominator = invert(T(degree, targetEpsilon.precision()));
        auto fDomain = Bounds<T>(LowerBound<T>(0, targetEpsilon.precision()), UpperBound<T>(1, targetEpsilon.precision()));
        auto fRange = function(fDomain);
        for (size_t i = 1; i <= degree; ++i) {
            auto x = i * denominator;
            auto interval = Bounds<T>(LowerBound<T>((i - 1) * denominator), UpperBound<T>(x));

            auto originalBounds = function(interval);
            auto polynomialBounds = this->evaluate(interval);

            polynomialBounds = refinement(polynomialBounds, fRange);

            auto maxError = mag((originalBounds - polynomialBounds) * 5 / 4);
/* 
            std::cout << interval << std::endl;
            std::cout << originalBounds << std::endl;
            std::cout << polynomialBounds << std::endl;
            std::cout << maxError << std::endl
                      << originalBounds - polynomialBounds << std::endl
                      << degree << std::endl
                      << std::endl;
 */
            if ((maxError > targetEpsilon).repr() >= LogicalValue::LIKELY)
                return false;

            _errorBounds.emplace_back(maxError);
        }

        return true;
    }

    static Bounds<T> invert(const Bounds<T>& x) {
        return T(1, x.precision()) / x;
    }
    static Bounds<T> invert(const T& x) {
        return T(1, x.precision()) / x;
    }

    std::vector<PositiveUpperBound<T>> _errorBounds{};
};

} // namespace Ariadne

#endif
