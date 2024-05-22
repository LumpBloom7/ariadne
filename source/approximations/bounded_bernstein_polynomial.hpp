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
class BoundedBernsteinPolynomial_impl : virtual public IBernsteinPolynomial<T> {
    using PR = T::PrecisionType;

  public:
    BoundedBernsteinPolynomial_impl(const IBernsteinPolynomialPtr<T> &polynomial, const std::function<Bounds<T>(Bounds<T>)> &function)
        : polynomial(polynomial) {
        computeErrorBounds(function, PositiveUpperBound<T>(T::inf(precision())));
    }
    BoundedBernsteinPolynomial_impl(const IBernsteinPolynomialPtr<T> &polynomial, const std::function<Bounds<T>(Bounds<T>)> &function, PositiveUpperBound<T> epsilon)
        : polynomial(polynomial) {
        computeErrorBounds(function, epsilon);
    }

    virtual Bounds<T> evaluate(const Bounds<T> &x) const override {
        return polynomial.evaluate(x);
    }

    virtual Bounds<T> evaluateDerivative(const Bounds<T> &x) const override {
        return polynomial.evaluateDerivative(x);
    }

    virtual DegreeType degree() const override {
        return polynomial.degree();
    };

    virtual PR precision() const override {
        return polynomial.precision();
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

    virtual std::shared_ptr<IBernsteinPolynomial<T>> asSharedPtr() const override {
        return std::make_shared<BoundedBernsteinPolynomial_impl>(*this);
    };

  protected:
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
    IBernsteinPolynomialPtr<T> polynomial;
};

template<typename T>
class BoundedBernsteinPolynomial : public IBernsteinPolynomialPtr<T> {
  public:
    using RND = T::RoundingModeType;
    using PR = T::PrecisionType;

    template<typename... TArgs>
    BoundedBernsteinPolynomial(TArgs... args) {
        (this->_ptr) = _ptr2 = std::make_shared<BoundedBernsteinPolynomial_impl<T>>(args...);
    }

    PositiveUpperBound<T> maximumError() const {
        return _ptr2->maximumError();
    }

    PositiveUpperBound<T> maximumErrorAt(const Bounds<T> &x) const {
        return _ptr2->maximumErrorAt(x);
    }

  protected:

  private:
    std::shared_ptr<BoundedBernsteinPolynomial_impl<T>> _ptr2{nullptr};
};

} // namespace Ariadne

#endif
