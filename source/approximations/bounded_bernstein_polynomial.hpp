#ifndef ARIADNE_BOUNDED_BERNSTEIN_POLYNOMIAL_HPP
#define ARIADNE_BOUNDED_BERNSTEIN_POLYNOMIAL_HPP

#include <cmath>
#include <functional>
#include <iostream>

#include "approximations/polynomial_approximation_interface.hpp"
#include "geometry/interval.hpp"

namespace Ariadne {
template<typename T> class BoundedPolynomialApproximation : virtual public IPolynomialApproximation<T> {
    using PR = T::PrecisionType;

    using TBounds = Bounds<T>;
    using TInterval = Interval<T>;

  public:
    template<typename Y, typename... TArgs>
    BoundedPolynomialApproximation(const std::function<TBounds(TBounds)> &function, TArgs... args)
        : originalPoly(std::make_shared<Y>(args...)) {
        computeErrorBounds(function);
    }
    template<typename Y, typename... TArgs>
    BoundedPolynomialApproximation(const std::function<TBounds(TBounds)> &function, PositiveUpperBound<T> epsilon,
                                   TArgs... args)
        : originalPoly(std::make_shared<Y>(args...)) {
        computeErrorBounds(function, epsilon);
    }

    BoundedPolynomialApproximation(const std::function<TBounds(TBounds)> &function,
                                   const std::shared_ptr<IPolynomialApproximation<T>> &approximationPtr,
                                   int subintervals = 1)
        : originalPoly(approximationPtr) {
        computeErrorBounds(function, subintervals);
    }

    BoundedPolynomialApproximation(const std::function<TBounds(TBounds)> &function,
                                   const std::shared_ptr<IPolynomialApproximation<T>> &approximationPtr,
                                   PositiveUpperBound<T> epsilon, int depth = 1, int subintervals = 1)
        : originalPoly(approximationPtr) {
        computeErrorBounds(function, epsilon, depth, subintervals);
    }

    BoundedPolynomialApproximation(const std::function<TBounds(TBounds)> &function,
                                   const std::shared_ptr<IPolynomialApproximation<T>> &approximationPtr,
                                   PositiveUpperBound<T> epsilon, const T &minimumIntervalWidth, int subintervals = 1)
        : originalPoly(approximationPtr) {
        computeErrorBounds(function, epsilon, minimumIntervalWidth, subintervals);
    }

    virtual TBounds range(const TInterval &x, Nat subInterval = 1) const override {
        return originalPoly->range(x, subInterval);
    }

    virtual TBounds evaluate(const TBounds &x) const override { return originalPoly->evaluate(x); }

    virtual TBounds evaluateDerivative(const TBounds &x) const override { return originalPoly->evaluateDerivative(x); }

    virtual DegreeType degree() const override { return originalPoly->degree(); }

    virtual PR precision() const override { return originalPoly->precision(); }

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

#ifndef _OPENMP
    void computeErrorBounds(const std::function<TBounds(TBounds)> &function, const PositiveUpperBound<T> &targetEpsilon,
                            int depth = 1, int subintervals = 1) {
        auto degree = this->degree();

        _errorBounds.clear();
        _errorBounds.reserve(degree);

        auto zero = TBounds(this->precision());
        auto degreeReciprocal = rec(TBounds(degree, zero.precision()));

        auto x = TBounds(zero.lower_raw(), degreeReciprocal.upper_raw());

        for (int i = 1; i <= degree; ++i) {
            auto error = computeErrorBounds(function, x, targetEpsilon, depth, subintervals);
            _errorBounds.emplace_back(error);

            if (decide(error > targetEpsilon))
                return;

            x += degreeReciprocal;
        }
    }

    PositiveUpperBound<T> computeErrorBounds(const std::function<TBounds(TBounds)> &function, const TBounds &interval,
                                             const PositiveUpperBound<T> &targetEpsilon, int depth = 1,
                                             int subintervals = 1) {
        auto actual = function(interval);
        auto error = mag(actual.upper_raw() - actual.lower_raw());

        if (decide(error <= targetEpsilon)) {
            auto predicted = evaluate(interval);
            error = mag(actual - predicted);
        }

        if (depth == 0 || subintervals == 1) // We're not making subintervals smaller, no point continuing
            return error;

        if (decide(error > targetEpsilon)) {
            auto intervalWidth = interval.upper_raw() - interval.lower_raw();

            auto step = intervalWidth / subintervals;
            auto subinterval = TBounds(interval.lower_raw(), UpperBound<T>(interval.lower_raw() + step));
            PositiveUpperBound<T> maxSubError = PositiveUpperBound<T>(this->precision());

            for (int i = 0; i < subintervals; ++i) {
                auto subError = computeErrorBounds(function, subinterval, targetEpsilon, depth - 1, subintervals);
                maxSubError = max(subError, maxSubError);
                subinterval += step;
            }

            error = maxSubError;
        }

        return error;
    }
#else
    // Will serve as the multithreaded cancellation token, but is not thread safe
    // itself
    bool threadCancellation = false;

    void computeErrorBounds(const std::function<TBounds(TBounds)> &function, const PositiveUpperBound<T> &targetEpsilon,
                            int depth = 1, int subintervals = 1) {
        threadCancellation = false;
        auto degree = this->degree();

        _errorBounds.clear();
        _errorBounds.reserve(degree);
        auto zero = TBounds(this->precision());

        for (int i = 0; i < degree; ++i)
            _errorBounds.emplace_back(zero);

        auto degreeReciprocal = rec(TBounds(degree, zero.precision()));
        auto x = TBounds(zero.lower_raw(), degreeReciprocal.upper_raw());

#pragma omp parallel for
        for (int i = 0; i < degree; ++i) {
            if (threadCancellation)
                continue;

            auto subInterval = x + (i * degreeReciprocal);
            auto error = computeErrorBounds(function, subInterval, targetEpsilon, depth, subintervals);
            _errorBounds[i] = error;

            if (decide(error > targetEpsilon))
                threadCancellation = true;
        }
    }

    PositiveUpperBound<T> computeErrorBounds(const std::function<TBounds(TBounds)> &function, const TBounds &interval,
                                             const PositiveUpperBound<T> &targetEpsilon, int depth = 1,
                                             int subintervals = 1) {
        if (threadCancellation)
            return PositiveUpperBound<T>(this->precision());

        auto actual = function(interval);
        auto error = mag(actual.upper_raw() - actual.lower_raw());

        if (decide(error <= targetEpsilon)) {
            auto predicted = evaluate(interval);
            error = mag(actual - predicted);
        }

        if (depth == 0 || subintervals == 1) // We're not making subintervals smaller, no point continuing
            return error;

        if (decide(error > targetEpsilon)) {
            auto intervalWidth = interval.upper_raw() - interval.lower_raw();

            auto step = intervalWidth / subintervals;
            auto firstSubInterval = TBounds(interval.lower_raw(), UpperBound<T>(interval.lower_raw() + step));

            std::vector<PositiveUpperBound<T>> subErrors = {};
            PositiveUpperBound<T> maxSubError = PositiveUpperBound<T>(this->precision());

            for (int i = 0; i < subintervals; ++i) {
                subErrors.emplace_back(maxSubError);
            }

#pragma omp parallel for
            for (int i = 0; i < subintervals; ++i) {
                auto subinterval = firstSubInterval + (i * step);
                auto subError = computeErrorBounds(function, subinterval, targetEpsilon, depth - 1, subintervals);

                subErrors[i] = subError;
            }

            for (int i = 0; i < subintervals; ++i) {
                maxSubError = max(subErrors[i], maxSubError);
            }

            error = maxSubError;
        }

        return error;
    }

#endif

#ifndef _OPENMP
    void computeErrorBounds(const std::function<TBounds(TBounds)> &function, int subintervals = 1) {
        auto degree = this->degree();

        _errorBounds.clear();
        _errorBounds.reserve(degree);

        auto zero = TBounds(this->precision());
        auto degreeReciprocal = rec(TBounds(degree, zero.precision()));

        auto x = TBounds(zero.lower_raw(), degreeReciprocal.upper_raw());

        for (int i = 1; i <= degree; ++i) {
            auto actual = function(x);
            auto predicted = this->range(toInterval(x), subintervals);

            auto errorUpperBound = mag(actual - predicted);

            _errorBounds.emplace_back(errorUpperBound);

            x += degreeReciprocal;
        }
    }
#else
    void computeErrorBounds(const std::function<TBounds(TBounds)> &function, int subintervals = 1) {
        auto degree = this->degree();

        _errorBounds.clear();
        _errorBounds.reserve(degree);

        auto zero = TBounds(this->precision());
        auto degreeReciprocal = rec(TBounds(degree, zero.precision()));

        for (int i = 0; i < degree; ++i)
            _errorBounds.emplace_back(zero);

        auto x = TBounds(zero.lower_raw(), degreeReciprocal.upper_raw());
#pragma omp parallel for
        for (int i = 0; i < degree; ++i) {
            auto interval = x + degreeReciprocal * i;
            auto actual = function(interval);
            auto predicted = this->range(this->toInterval(interval), subintervals);

            auto errorUpperBound = mag(actual - predicted);

            _errorBounds[i] = errorUpperBound;
        }
    }
#endif

  private:
    std::vector<PositiveUpperBound<T>> _errorBounds{};
    std::shared_ptr<IPolynomialApproximation<T>> originalPoly;
};

} // namespace Ariadne

#endif
