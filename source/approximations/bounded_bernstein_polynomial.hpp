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
class BoundedPolynomialApproximation : public IPolynomialApproximation<T> {
    using PR = T::PrecisionType;

    using TBounds = Bounds<T>;

  public:
    template<typename Y, typename... TArgs>
    BoundedPolynomialApproximation(const std::function<TBounds(TBounds)> &function, TArgs... args)
        : originalPoly(std::make_shared<Y>(args...)) {
        computeErrorBounds(function);
    }
    template<typename Y, typename... TArgs>
    BoundedPolynomialApproximation(const std::function<TBounds(TBounds)> &function, PositiveUpperBound<T> epsilon, TArgs... args)
        : originalPoly(std::make_shared<Y>(args...)) {
        computeErrorBounds(function, epsilon);
    }

    BoundedPolynomialApproximation(const std::function<TBounds(TBounds)> &function, const std::shared_ptr<IPolynomialApproximation<T>> &approximationPtr, int subintervals = 1)
        : originalPoly(approximationPtr) {
        computeErrorBounds(function, subintervals);
    }

    BoundedPolynomialApproximation(const std::function<TBounds(TBounds)> &function, const std::shared_ptr<IPolynomialApproximation<T>> &approximationPtr, PositiveUpperBound<T> epsilon, int depth = 1, int subintervals = 1)
        : originalPoly(approximationPtr) {
        computeErrorBounds(function, epsilon, depth, subintervals);
    }

    virtual TBounds evaluate(const TBounds &x, int subInterval = 1) const override {
        return originalPoly->evaluate(x, subInterval);
    }
    virtual TBounds evaluateRaw(const TBounds &x) const override {
        return originalPoly->evaluateRaw(x);
    }

    virtual TBounds evaluateDerivative(const TBounds &x) const override {
        return originalPoly->evaluateDerivative(x);
    }

    virtual DegreeType degree() const override {
        return originalPoly->degree();
    }

    virtual PR precision() const override {
        return originalPoly->precision();
    }

    PositiveUpperBound<T> maximumError() const {
        auto maximum = _errorBounds[0];

        for (int i = 1; i < _errorBounds.size(); ++i) {
            maximum = max(maximum, _errorBounds[i]);

            if (is_inf(_errorBounds[i].raw())) {
                std::cout << i << std::endl;
            }
        }

        return maximum;
    }

    PositiveUpperBound<T> maximumErrorAt(const TBounds &x) const {
        auto degree = this->degree();
        auto maximum = PositiveUpperBound<T>(x.precision());

        auto xr = TBounds(x.precision());
        auto degreeReciprocal = rec(TBounds(degree, x.precision()));

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

    void computeErrorBounds(const std::function<TBounds(TBounds)> &function, const PositiveUpperBound<T> &targetEpsilon, int depth = 1, int subintervals = 1) {
        auto degree = this->degree();

        _errorBounds.clear();
        _errorBounds.reserve(degree);

        auto zero = TBounds(this->precision());
        auto degreeReciprocal = rec(TBounds(degree, zero.precision()));

        auto x = TBounds(zero.lower_raw(), degreeReciprocal.upper_raw());

        for (int i = 1; i <= degree; ++i) {
            auto error = computeErrorBounds(function, x, targetEpsilon, depth - 1, subintervals);
            _errorBounds.emplace_back(error);

            if (decide(error > targetEpsilon))
                return;

            x += degreeReciprocal;
        }
    }

    PositiveUpperBound<T> computeErrorBounds(const std::function<TBounds(TBounds)> &function, const TBounds &interval, const PositiveUpperBound<T> &targetEpsilon, int depth = 1, int subintervals = 1) {
        if (depth == 0) { // Base case
            auto actual = function(interval);
            auto predicted = this->evaluate(interval, subintervals);

            auto error = mag(actual - predicted);

            return error;
        }

        auto step = (interval.upper_raw() - interval.lower_raw()) / subintervals;
        auto x = TBounds(interval.lower_raw(), UpperBound<T>(interval.lower_raw() + step));
        PositiveUpperBound<T> maxError = PositiveUpperBound<T>(this->precision());

        for (int i = 0; i < subintervals; ++i) {
            auto actual = function(x);
            auto error = mag(actual.upper_raw() - actual.lower_raw());

            if (decide(error <= targetEpsilon)) {
                auto predicted = this->evaluate(x);
                error = mag(actual - predicted);
            }

            if (decide(error > targetEpsilon)) {
                error = computeErrorBounds(function, x, targetEpsilon, depth - 1, subintervals);
            }

            maxError = max(error, maxError);

            if (decide(error > targetEpsilon)) {
                break;
            }
            x += step;
        }

        return maxError;
    }

    void computeErrorBounds(const std::function<TBounds(TBounds)> &function, int subintervals = 1) {
        auto degree = this->degree();

        _errorBounds.clear();
        _errorBounds.reserve(degree);

        auto zero = TBounds(this->precision());
        auto degreeReciprocal = rec(TBounds(degree, zero.precision()));

        auto x = TBounds(zero.lower_raw(), degreeReciprocal.upper_raw());

        for (int i = 1; i <= degree; ++i) {
            auto actual = function(x);
            auto predicted = this->evaluate(x, subintervals);

            auto errorUpperBound = mag(actual - predicted);

            _errorBounds.emplace_back(errorUpperBound);

            x += degreeReciprocal;
        }
    }

  private:
    std::vector<PositiveUpperBound<T>> _errorBounds{};
    std::shared_ptr<IPolynomialApproximation<T>> originalPoly;
};

} // namespace Ariadne

#endif
