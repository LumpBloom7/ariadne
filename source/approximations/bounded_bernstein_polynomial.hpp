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
  public:
    BoundedBernsteinPolynomial(const std::function<Bounds<T>(Bounds<T>)> &function, int degree, size_t maxPrecision)
        : BernsteinPolynomial<T>(function, degree, maxPrecision) {
        computeErrorBounds(function, PositiveUpperBound<T>(T::inf(MultiplePrecision(maxPrecision))));
    }

    BoundedBernsteinPolynomial(const std::function<Bounds<T>(Bounds<T>)> &function, PositiveUpperBound<T> targetEpsilon)
        : BernsteinPolynomial<T>({}) {
        int degree = 1;
        do {
            this->generateCoefficients(function, degree, targetEpsilon.precision());
            degree *= 2;
        } while (!computeErrorBounds(function, targetEpsilon));
    }

    PositiveUpperBound<T> maximumErrorAt(Bounds<T> x) {
        auto degree = (this->_coefficients).size() - 1;
        Bounds<T> denominator = Bounds<T>(T(Approximation<T>(1 / static_cast<float>(degree), x.precision())));
        auto maximum = PositiveUpperBound<T>(x.precision());

        for (size_t i = 0; i <= degree; ++i) {
            auto xl = i * denominator;
            auto xr = (i + 1) * denominator;

            if ((xl > x.upper() || xr < x.lower()).repr() >= LogicalValue::LIKELY)
                continue;

            std::cout << "t" << std::endl;

            maximum = max(maximum, _errorBounds[i]);
        }
        std::cout << std::endl;

        return maximum;
    }

  private:
    bool computeErrorBounds(const std::function<Bounds<T>(Bounds<T>)> &function, PositiveUpperBound<T> targetEpsilon) {
        auto degree = (this->_coefficients).size() - 1;
        _errorBounds = {};
        T denominator = T(Approximation<T>(1 / static_cast<float>(degree), targetEpsilon.precision()));

        for (size_t i = 1; i <= degree; ++i) {
            auto x = i * denominator;
            auto interval = Bounds<T>(LowerBound<T>((i - 1) * denominator), UpperBound<T>(x));

            auto originalBounds = function(interval);
            auto polynomialBounds = this->DeCasteljau(interval);

            auto maxError = mag(originalBounds - polynomialBounds);

            std::cout << interval << std::endl;
            std::cout << originalBounds << std::endl;
            std::cout << polynomialBounds << std::endl;
            std::cout << maxError << std::endl
                      << std::endl;

            if ((maxError > targetEpsilon).repr() >= LogicalValue::LIKELY)
                return false;

            _errorBounds.emplace_back(maxError);
        }

        return true;
    }

    std::vector<PositiveUpperBound<T>> _errorBounds;
};

} // namespace Ariadne

#endif
