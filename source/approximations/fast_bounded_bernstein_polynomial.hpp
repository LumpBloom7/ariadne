#ifndef ARIADNE_FAST_BOUNDED_BERNSTEIN_POLYNOMIAL_HPP
#define ARIADNE_FAST_BOUNDED_BERNSTEIN_POLYNOMIAL_HPP

#include <cmath>
#include <functional>
#include <iostream>

#include "approximations/polynomial_approximation_interface.hpp"
#include "geometry/interval.hpp"
#include "numeric/float_bounds.hpp"

namespace Ariadne {
template<typename T> class FastBoundedPolynomialApproximation : virtual public IPolynomialApproximation<T> {
    using PR = T::PrecisionType;

    using TBounds = Bounds<T>;
    using TInterval = Interval<T>;

  public:
    template<typename Y, typename... TArgs>
    FastBoundedPolynomialApproximation(const std::function<TBounds(TBounds)> &function, TArgs... args)
        : originalPoly(std::make_shared<Y>(args...)) {
        computeErrorBounds(function);
    }
    template<typename Y, typename... TArgs>
    FastBoundedPolynomialApproximation(const std::function<TBounds(TBounds)> &function, PositiveUpperBound<T> epsilon,
                                       TArgs... args)
        : originalPoly(std::make_shared<Y>(args...)) {
        computeErrorBounds(function, epsilon);
    }

    FastBoundedPolynomialApproximation(const std::function<TBounds(TBounds)> &function,
                                       const std::shared_ptr<IPolynomialApproximation<T>> &approximationPtr,
                                       int secantInterations = 1)
        : originalPoly(approximationPtr) {
        computeErrorBounds(function, secantInterations);
    }

    FastBoundedPolynomialApproximation(const std::function<TBounds(TBounds)> &function,
                                       const std::shared_ptr<IPolynomialApproximation<T>> &approximationPtr,
                                       PositiveUpperBound<T> epsilon, int secantInterations = 1)
        : originalPoly(approximationPtr) {
        computeErrorBounds(function, epsilon, secantInterations);
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

    void computeErrorBounds(const std::function<TBounds(TBounds)> &function, int secantIterations = 1) {
        auto degree = this->degree();

        _errorBounds.clear();
        _errorBounds.reserve(degree);

        auto zero = TBounds(this->precision());
        auto degreeReciprocal = rec(TBounds(degree, zero.precision()));

        auto x = TBounds(zero.lower_raw(), degreeReciprocal.upper_raw());

        for (int i = 1; i <= degree; ++i) {
            auto actual = function(x);

            auto critX = secantMethod(x, secantIterations);

            auto left = evaluate(x.lower_raw());
            auto right = evaluate(x.upper_raw());

            auto predicted = unionOf(left, right);

            if (models(x, critX.value_raw())) {
                auto crit = evaluate(critX);
                predicted = unionOf(predicted, crit);
            }

            auto errorUpperBound = mag(actual - predicted);

            _errorBounds.emplace_back(errorUpperBound);

            x += degreeReciprocal;
        }
    }

    void computeErrorBounds(const std::function<TBounds(TBounds)> &function, PositiveUpperBound<T> targetEpsilon,
                            int secantIterations = 1) {
        auto degree = this->degree();

        _errorBounds.clear();
        _errorBounds.reserve(degree);

        auto zero = TBounds(this->precision());
        auto degreeReciprocal = rec(TBounds(degree, zero.precision()));

        auto x = TBounds(zero.lower_raw(), degreeReciprocal.upper_raw());

        for (int i = 1; i <= degree; ++i) {
            auto actual = function(x);

            auto critX = secantMethod(x, secantIterations);

            auto left = evaluate(x.lower_raw());
            auto right = evaluate(x.upper_raw());

            auto predicted = unionOf(left, right);

            if (models(x, critX.value_raw())) {
                auto crit = evaluate(critX);
                predicted = unionOf(predicted, crit);
            }

            auto errorUpperBound = mag(actual - predicted);
            _errorBounds.emplace_back(errorUpperBound);

            if (decide(errorUpperBound > targetEpsilon)) {
                return;
            }

            x += degreeReciprocal;
        }
    }

  private:
    const Bounds<T> secantMethod(const Bounds<T> &x, const PositiveUpperBound<T> &targetEpsilon) const {
        T s[]{x.lower_raw(), x.upper_raw()};

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
        T s[]{x.lower_raw(), x.upper_raw()};

        auto leftval = this->evaluateDerivative(s[0]);

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

        auto rightDeriv = this->evaluateDerivative(right);

        auto res = (right - rightDeriv * ((right - left) / (rightDeriv - leftDeriv))).value();

        if (is_nan(res))
            return false;

        x[0] = x[1];
        x[1] = res;
        leftDeriv = rightDeriv;
        return true;
    }

    Bounds<T> unionOf(Bounds<T> first, Bounds<T> second) {
        return Bounds<T>(min(first.lower_raw(), second.lower_raw()), max(first.upper_raw(), second.upper_raw()));
    }

    std::vector<PositiveUpperBound<T>> _errorBounds{};
    std::shared_ptr<IPolynomialApproximation<T>> originalPoly;
};

} // namespace Ariadne

#endif
