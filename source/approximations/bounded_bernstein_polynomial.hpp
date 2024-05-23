#ifndef ARIADNE_BOUNDED_BERNSTEIN_POLYNOMIAL_HPP
#define ARIADNE_BOUNDED_BERNSTEIN_POLYNOMIAL_HPP

#include <cmath>
#include <deque>
#include <exception>
#include <functional>
#include <iostream>
#include <optional>
#include <stack>

#include "approximations/bernstein_polynomial.hpp"
namespace Ariadne {
template<typename T>
class BoundedBernsteinPolynomial : public BernsteinPolynomial<T> {
    using PR = T::PrecisionType;

  public:
    BoundedBernsteinPolynomial(const std::function<Bounds<T>(Bounds<T>)> &function, const BernsteinPolynomial<T> &original) : BernsteinPolynomial<T>(original) {
        computeErrorBounds(function, PositiveUpperBound<T>(T::inf(this->precision())));
    }
    BoundedBernsteinPolynomial(const std::function<Bounds<T>(Bounds<T>)> &function, const BernsteinPolynomial<T> &original, PositiveUpperBound<T> epsilon) : BernsteinPolynomial<T>(original) {
        computeErrorBounds(function, epsilon);
    }

    BoundedBernsteinPolynomial(const std::function<Bounds<T>(Bounds<T>)> &function, DegreeType degree, PR precision, int secantIters = 5) : BoundedBernsteinPolynomial(function, PositiveUpperBound<T>(T::inf(precision)), degree, precision, secantIters) {
    }
    BoundedBernsteinPolynomial(const std::function<Bounds<T>(Bounds<T>)> &function, PositiveUpperBound<T> epsilon, DegreeType degree, PR precision, int secantIters = 5) : BernsteinPolynomial<T>(function, degree, precision, secantIters) {
        computeErrorBounds(function, epsilon);
    }

    PositiveUpperBound<T> maximumError() const {
        auto maximum = _errorBounds[0];

        for (int i = 1; i < _errorBounds.size(); ++i) {
            maximum = max(maximum, _errorBounds[i]);
        }

        return maximum;
    }

    PositiveUpperBound<T> maximumErrorAt(const Bounds<T> &x) const {
        auto degree = this->degree();
        auto maximum = PositiveUpperBound<T>(x.precision());

        auto xr = Bounds<T>(x.precision());
        auto degreeReciprocal = rec(Bounds<T>(degree, x.precision()));

        for (size_t i = 0; i <= degree; ++i) {
            xr += degreeReciprocal;

            if ((x.lower_raw() >= xr).repr() >= LogicalValue::LIKELY)
                continue;

            maximum = max(maximum, _errorBounds[i]);

            if ((x.upper_raw() <= xr).repr() >= LogicalValue::LIKELY)
                break;
        }

        return maximum;
    }

    void computeErrorBounds(const std::function<Bounds<T>(Bounds<T>)> &function, const PositiveUpperBound<T> &targetEpsilon) {
        auto degree = this->degree();

        _errorBounds.clear();
        _errorBounds.reserve(degree);

        auto zero = Bounds<T>(targetEpsilon.precision());
        auto degreeReciprocal = rec(Bounds<T>(degree, zero.precision()));

        auto x = Bounds<T>(zero.lower_raw(), degreeReciprocal.upper_raw());

        for (int i = 1; i <= degree; ++i) {
            auto actual = function(x);
            auto predicted = this->evaluate(x);

            auto errorUpperBound = mag(actual - predicted);

            _errorBounds.emplace_back(errorUpperBound);

            if ((errorUpperBound > targetEpsilon).repr() >= LogicalValue::INDETERMINATE)
                return;

            x += degreeReciprocal;
        }
    }
  private:
    std::vector<PositiveUpperBound<T>> _errorBounds{};
};

} // namespace Ariadne

#endif
