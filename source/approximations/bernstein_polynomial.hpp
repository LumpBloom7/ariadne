#ifndef ARIADNE_BERNSTEIN_POLYNOMIAL_HPP
#define ARIADNE_BERNSTEIN_POLYNOMIAL_HPP

#include <cmath>
#include <functional>
#include <iostream>

#include "numeric/bits.hpp"
#include "numeric/dyadic.hpp"
#include "numeric/float_bounds.hpp"
#include "numeric/floatmp.hpp"
#include "numeric/int.hpp"
#include "numeric/rational.hpp"
#include "numeric/real.hpp"
#include "numeric/rounding.hpp"
#include "numeric/validated_real.hpp"
#include "utility/standard.hpp"

namespace Ariadne {

class BernsteinPolynomial {
  public:
    BernsteinPolynomial(const std::function<double(double)>& function, int degree);
    FloatMPBounds evaluate(double x, size_t effort);

  private:
    FloatMPBounds binomialCoefficient(const int n, const int k);
    std::vector<ExactDouble> computeCoefficients();
    FloatMPBounds bernsteinBasisPolynomialFor(int v, double x, size_t effort);

    static FloatMPBounds factorial(int x);

    const std::function<double(double)>& _function;
    const int _degree;
    std::vector<ExactDouble> _coefficients;

    // Caches to reduce computation time spent on things not dependent on x
    // Using these allows the computation time to be O(1) amortised
    static std::vector<FloatMPBounds> _factorialCache;
    static double _mantissaBitsCount;
    std::vector<FloatMPBounds> _binomialCache;

    static FloatMPBounds toPrecision(const FloatMPBounds& x, const MultiplePrecision& pr)
    {
        return FloatMPBounds(x, pr);
    }
};

} // namespace Ariadne

#endif // ARIADNE_BERNSTEIN_POLYNOMIAL_HPP
