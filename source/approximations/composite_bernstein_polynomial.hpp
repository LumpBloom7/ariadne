#pragma once
#include "bernstein_polynomial.hpp"

namespace Ariadne {

template<typename T>
class CompositeBernsteinPolynomial : public IBernsteinPolynomial<T> {
    using BType = IBernsteinPolynomial<T>;
    using PR = T::PrecisionType;

  public:
    CompositeBernsteinPolynomial(std::initializer_list<BType> list) {
        for (BType& p: list)
            polynomials.emplace_back(p);
    };

    virtual Bounds<T> evaluate(const Bounds<T> &x) const override {
        Bounds<T> sum = Bounds<T>(x.precision());

        for (BType p: polynomials)
            sum += p.evaluate(x);

        return sum;
    }

    virtual Bounds<T> evaluateDerivative(const Bounds<T> &x) const override {
        Bounds<T> sum = Bounds<T>(x.precision());

        for (BType& p: polynomials)
            sum += p.evaluateDerivative(x);

        return sum;
    }

    virtual DegreeType degree() const override {
        DegreeType maxDegree = 0;

        for (BType& p: polynomials) {
            if (p.degree() > maxDegree)
                maxDegree = p.degree();
        }

        return maxDegree;
    }

    virtual PR precision() const override {
        if (polynomials.size() == 0)
            return PR(0);

        PR minPrecision = polynomials[0].precision();

        for (BType& p: polynomials)
            minPrecision = min(p.precision(), minPrecision);

        return minPrecision;
    }

    void add(const BType &element) {
        polynomials.emplace_back(element);
    }

    CompositeBernsteinPolynomial<T> &operator+=(const BType &p) {
        add(p);
        return *this;
    }
    CompositeBernsteinPolynomial<T> &operator-=(const BType &p) {
        add(BernsteinPolynmomialNegation<T>(p));
        return *this;
    }

  private:
    std::vector<BType> polynomials = {};
};

template<typename T>
CompositeBernsteinPolynomial<T> operator+(const IBernsteinPolynomial<T> &p1, const IBernsteinPolynomial<T> &p2) { return CompositeBernsteinPolynomial({p1, p2}); }

template<typename T>
CompositeBernsteinPolynomial<T> operator-(const IBernsteinPolynomial<T> &p1, const IBernsteinPolynomial<T> &p2) { return CompositeBernsteinPolynomial({p1, BernsteinPolynmomialNegation<T>(p2)}); }
} // namespace Ariadne
