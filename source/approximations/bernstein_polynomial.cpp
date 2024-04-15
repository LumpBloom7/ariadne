#include "approximations/bernstein_polynomial.hpp"

#include "numeric/rational.hpp"

namespace Ariadne {

BernsteinPolynomial::BernsteinPolynomial(const std::function<double(double)>& function, int degree) : _function{function}, _degree{degree} {
    _coefficients = computeCoefficients();
    _binomialCache = {};
}

Integer BernsteinPolynomial::binomialCoefficient(const int n, const int k) {
    // Like the factorial function, this one uses a cache as well. The precision should always be larger-than or equal to what is actually needed.
    // This is stored not on a global level though.

    if (_binomialCache.size() > k)
        return _binomialCache[k];

    // FloatMP preserves precision of the least precise factor during multiplication of two FloatMPs
    // We try to resolve this by first "casting" them to a higher precision level (which is overestimated)
    auto a = Factorials::get(n);
    auto b = Factorials::get(k), c = Factorials::get(n - k);

    auto res = (a / (b * c)).get_num(); // This will always be an integer

    return _binomialCache.emplace_back(res);
}

std::vector<ExactDouble> BernsteinPolynomial::computeCoefficients() {
    std::vector<ExactDouble> coefficients = std::vector<ExactDouble>(_degree + 1);

    for (size_t i = 0; i <= _degree; ++i)
        coefficients[i] = ExactDouble(_function(i / static_cast<float>(_degree)));

    return coefficients;
}

FloatMPBounds BernsteinPolynomial::bernsteinBasisPolynomialFor(int v, double x, size_t effort) {
    auto X = FloatMPBounds(exact(x), MultiplePrecision(effort));

    return binomialCoefficient(_degree, v) * pow(X, v) * pow(1 - X, _degree - v);
}

FloatMPBounds BernsteinPolynomial::evaluate(double x, size_t effort) {
    FloatMPBounds sum = FloatMPBounds(exact(0), MultiplePrecision(effort));

    if (x < 0 || x > 1)
        return sum;

    for (size_t i = 0; i <= _degree; ++i) {
        auto bp = bernsteinBasisPolynomialFor(i, x, effort);
        sum += _coefficients[i] * bp;
    }

    return sum;
}
} // namespace Ariadne
