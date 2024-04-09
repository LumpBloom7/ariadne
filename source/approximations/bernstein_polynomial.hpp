#ifndef ARIADNE_BERNSTEIN_POLYNOMIAL_HPP
#define ARIADNE_BERNSTEIN_POLYNOMIAL_HPP

#include "utility/standard.hpp"
#include <functional>
#include <cmath>
#include <iostream>
#include "numeric/int.hpp"

#include "numeric/floatmp.hpp"
#include "numeric/bits.hpp"
#include "numeric/float_bounds.hpp"
#include "numeric/rounding.hpp"

namespace Ariadne {

class BernsteinPolynomial {
public:
    BernsteinPolynomial(const std::function<double(double)>& function, size_t degree);
    Bounds<FloatMP> evaluate(double x);

private:
    Bounds<FloatMP> binomialCoefficient(const size_t n, const size_t k);
    std::vector<ExactDouble> computeCoefficients();
    Bounds<FloatMP> bernsteinBasisPolynomialFor(size_t v, double x);

    static Bounds<FloatMP> factorial(size_t x);

    const std::function<double(double)>& _function;
    const size_t _degree;
    std::vector<ExactDouble> _coefficients;

    // Caches to reduce computation time spent on things not dependent on x
    // Using these allows the computation time to be O(1) amortised
    static std::vector<Bounds<FloatMP>> _factorialCache;
    std::vector<Bounds<FloatMP>> _binomialCache;
};

}

#endif // ARIADNE_BERNSTEIN_POLYNOMIAL_HPP
