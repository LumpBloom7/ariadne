#include "approximations/bernstein_polynomial.hpp"

#include <iostream>

#include "numeric/rational.hpp"
#include "utility/hash_numeric.hpp"
#include "utility/hash_tuple.hpp"

namespace Ariadne {
namespace {
class BinomialCoefficients : Cache<Integer, std::tuple<Nat64, Nat64>> {
  public:
    using Cache<Integer, std::tuple<Nat64, Nat64>>::get;
    Integer operator()(const int n, const int m) { return get(n, m); }
    Integer get(const int n, const int m) { return get(std::make_tuple(n, m)); }

  protected:
    Integer compute(std::tuple<Nat64, Nat64> args) {
        Nat64 n = std::get<0>(args);
        Nat64 k = std::get<1>(args);

        auto a = Factorials::get(n);
        auto b = Factorials::get(k), c = Factorials::get(n - k);

        auto res = (a / (b * c)).get_num(); // This will always be an integer

        std::cout << res << std::endl;
        return res;
    }
};

static BinomialCoefficients binomialCoefficients = BinomialCoefficients();
} // namespace

BernsteinPolynomial::BernsteinPolynomial(const std::function<double(double)>& function, int degree) : _function{function}, _degree{degree} {
    _coefficients = computeCoefficients();
}

std::vector<ExactDouble> BernsteinPolynomial::computeCoefficients() {
    std::vector<ExactDouble> coefficients = std::vector<ExactDouble>(_degree + 1);

    for (size_t i = 0; i <= _degree; ++i)
        coefficients[i] = ExactDouble(_function(i / static_cast<float>(_degree)));

    return coefficients;
}

FloatMPBounds BernsteinPolynomial::bernsteinBasisPolynomialFor(int v, double x, size_t effort) {
    auto X = FloatMPBounds(exact(x), MultiplePrecision(effort));

    return binomialCoefficients(_degree, v) * pow(X, v) * pow(1 - X, _degree - v);
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
