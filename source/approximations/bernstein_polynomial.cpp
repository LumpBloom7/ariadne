#include "approximations/bernstein_polynomial.hpp"

namespace Ariadne {

BernsteinPolynomial::BernsteinPolynomial(const std::function<double(double)>& function, int degree) : _function{function}, _degree{degree}
{
    _coefficients = computeCoefficients();
    _binomialCache = {};
}

FloatMPBounds BernsteinPolynomial::factorial(int x)
{
    // This factorial function uses a cache to avoid recomputation of factors.
    // *In theory* all generated factorials should be 100% correct as I make sure that more than enough bits are reserved for the mantissa.
    // The cache is global, so all instances of BernsteinPolynomials will refer to the same cache.

    if (_factorialCache.size() > x)
        return _factorialCache[x];

    auto res = _factorialCache.back();

    // Recursion causes a stackoverflow on higher x values
    // For-loop my beloved gets to shine again.
    for (int i = _factorialCache.size(); i <= x; ++i) {
        _mantissaBitsCount += log2(i);
        int numbits = static_cast<int>((_mantissaBitsCount) + 1);

        res = FloatMPBounds(res, MultiplePrecision(numbits)) * i;

        std::cout << res << res.precision() << std::endl;

        _factorialCache.emplace_back(res);
    }

    return res;
}

FloatMPBounds BernsteinPolynomial::binomialCoefficient(const int n, const int k)
{
    // Like the factorial function, this one uses a cache as well. The precision should always be larger-than or equal to what is actually needed.
    // This is stored not on a global level though.

    if (_binomialCache.size() > k)
        return _binomialCache[k];

    // FloatMP preserves precision of the least precise factor during multiplication of two FloatMPs
    // We try to resolve this by first "casting" them to a higher precision level (which is overestimated)
    auto a = factorial(n);
    auto b = FloatMPBounds(factorial(k), a.precision()), c = FloatMPBounds(factorial(n - k), a.precision());

    return _binomialCache.emplace_back(a / (b * c));
}

std::vector<ExactDouble> BernsteinPolynomial::computeCoefficients()
{
    std::vector<ExactDouble> coefficients = std::vector<ExactDouble>(_degree + 1);

    for (size_t i = 0; i <= _degree; ++i)
        coefficients[i] = ExactDouble(_function(i / static_cast<float>(_degree)));

    return coefficients;
}

FloatMPBounds BernsteinPolynomial::bernsteinBasisPolynomialFor(int v, double x, size_t effort)
{
    auto X = FloatMPBounds(exact(x), MultiplePrecision(effort));

    return binomialCoefficient(_degree, v) * pow(X, v) * pow(1 - X, _degree - v);
}

FloatMPBounds BernsteinPolynomial::evaluate(double x, size_t effort)
{
    FloatMPBounds sum = FloatMPBounds(exact(0), MultiplePrecision(effort));

    if (x < 0 || x > 1)
        return sum;

    for (size_t i = 0; i <= _degree; ++i) {
        auto bp = bernsteinBasisPolynomialFor(i, x, effort);
        sum += _coefficients[i] * bp;
    }

    return sum;
}

std::vector<FloatMPBounds> BernsteinPolynomial::_factorialCache = {FloatMPBounds(exact(1), MultiplePrecision(1)), FloatMPBounds(exact(1), MultiplePrecision(1))};
double BernsteinPolynomial::_mantissaBitsCount = log2(1);

} // namespace Ariadne
