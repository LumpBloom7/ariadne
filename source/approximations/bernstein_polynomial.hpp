#ifndef ARIADNE_BERNSTEIN_POLYNOMIAL_HPP
#define ARIADNE_BERNSTEIN_POLYNOMIAL_HPP

#include <cmath>
#include <functional>
#include <iostream>

#include "numeric/float_bounds.hpp"
#include "numeric/floatmp.hpp"
#include "numeric/integer.hpp"
#include "utility/factorials.hpp"
#include "utility/standard.hpp"
#include "utility/cache.hpp"

namespace Ariadne {

class BernsteinPolynomial {
  public:
    BernsteinPolynomial(const std::function<double(double)>& function, int degree);
    FloatMPBounds evaluate(double x, size_t effort);

  private:
    std::vector<ExactDouble> computeCoefficients();
    FloatMPBounds bernsteinBasisPolynomialFor(int v, double x, size_t effort);

    const std::function<double(double)>& _function;
    const int _degree;
    std::vector<ExactDouble> _coefficients;
};

} // namespace Ariadne

#endif // ARIADNE_BERNSTEIN_POLYNOMIAL_HPP
