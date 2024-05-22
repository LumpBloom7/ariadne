#ifndef ARIADNE_BERNSTEIN_POLYNOMIAL_SHARED_HPP
#define ARIADNE_BERNSTEIN_POLYNOMIAL_SHARED_HPP

#include "utility/cache.hpp"
#include "utility/factorials.hpp"

namespace Ariadne {
class IBernsteinPolynomialBase {
  protected:
    static SimpleCache<Integer, Nat64, Nat64> binomialCoefficients;
};

template<typename T>
class IBernsteinPolynomial {
  protected:
    using RND = T::RoundingModeType;
    using PR = T::PrecisionType;

  public:
    ~IBernsteinPolynomial() {}
    virtual Bounds<T> evaluate(const Bounds<T> &x) const = 0;
    virtual Bounds<T> evaluateDerivative(const Bounds<T> &x) const = 0;

    virtual DegreeType degree() const = 0;
    virtual PR precision() const = 0;

    Bounds<T> operator()(const Bounds<T> &x) const {
        return evaluate(x);
    }

    virtual std::shared_ptr<IBernsteinPolynomial<T>> asSharedPtr() const = 0;
};

template<typename T>
class IBernsteinPolynomialPtr : virtual public IBernsteinPolynomial<T> {
    using PR = T::PrecisionType;

  public:
    IBernsteinPolynomialPtr() = default;
    IBernsteinPolynomialPtr(IBernsteinPolynomial<T> p) : _ptr(p.asSharedPtr()) {}

    virtual Bounds<T> evaluate(const Bounds<T> &x) const override {
        return _ptr->evaluate(x);
    }

    virtual Bounds<T> evaluateDerivative(const Bounds<T> &x) const override {
        return _ptr->evaluateDerivative(x);
    }

    virtual DegreeType degree() const override {
        return _ptr->degree();
    }

    virtual PR precision() const override {
        return _ptr->precision();
    }

    virtual std::shared_ptr<IBernsteinPolynomial<T>> asSharedPtr() const override {
        return _ptr;
    };

  protected:
    std::shared_ptr<IBernsteinPolynomial<T>> _ptr{nullptr};
};
} // namespace Ariadne

#endif
