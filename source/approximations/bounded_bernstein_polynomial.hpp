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

  public:
    template<typename Y, typename... TArgs>
    BoundedPolynomialApproximation(const std::function<Bounds<T>(Bounds<T>)> &function, TArgs... args)
        : originalPoly(std::make_shared<Y>(args...)) {
        computeErrorBounds(function);
    }
    template<typename Y, typename... TArgs>
    BoundedPolynomialApproximation(const std::function<Bounds<T>(Bounds<T>)> &function, PositiveUpperBound<T> epsilon, TArgs... args)
        : originalPoly(std::make_shared<Y>(args...)) {
        computeErrorBounds(function, epsilon);
    }

    BoundedPolynomialApproximation(const std::function<Bounds<T>(Bounds<T>)> &function, const std::shared_ptr<IPolynomialApproximation<T>> &approximationPtr)
        : originalPoly(approximationPtr) {
        computeErrorBounds(function);
    }

    BoundedPolynomialApproximation(const std::function<Bounds<T>(Bounds<T>)> &function, const std::shared_ptr<IPolynomialApproximation<T>> &approximationPtr, PositiveUpperBound<T> epsilon)
        : originalPoly(approximationPtr) {
        computeErrorBounds(function, epsilon);
    }

    virtual Bounds<T> evaluate(const Bounds<T> &x) const override {
        return originalPoly->evaluate(x);
    }
    virtual Bounds<T> evaluateRaw(const Bounds<T> &x) const override {
        return originalPoly->evaluateRaw(x);
    }

    virtual Bounds<T> evaluateDerivative(const Bounds<T> &x) const override {
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
            auto predicted = this->evaluateRaw(x);

            auto errorUpperBound = mag(actual - predicted);

            _errorBounds.emplace_back(errorUpperBound);

            if ((errorUpperBound > targetEpsilon).repr() >= LogicalValue::INDETERMINATE)
                return;

            x += degreeReciprocal;
        }
    }

    void computeErrorBounds(const std::function<Bounds<T>(Bounds<T>)> &function) {
        auto degree = this->degree();

        _errorBounds.clear();
        _errorBounds.reserve(degree);

        auto zero = Bounds<T>(this->precision());
        auto degreeReciprocal = rec(Bounds<T>(degree, zero.precision()));

        auto x = Bounds<T>(zero.lower_raw(), degreeReciprocal.upper_raw());

        for (int i = 1; i <= degree; ++i) {
            auto actual = function(x);
            auto predicted = this->evaluateRaw(x);

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
