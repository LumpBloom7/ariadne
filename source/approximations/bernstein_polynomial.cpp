#include "approximations/bernstein_polynomial.hpp"

namespace Ariadne {

Bounds<FloatMP>  BernsteinPolynomial::factorial(size_t x) {
    if(_factorialCache.size() > x)
        return _factorialCache[x];

    Ariadne::Bounds<FloatMP> res = _factorialCache.back();

    for(size_t i = _factorialCache.size(); i <= x; ++i) {
        res = res * i;
        _factorialCache.emplace_back(res);
    }

    return res;
}

Bounds<FloatMP> BernsteinPolynomial::binomialCoefficient(const size_t n, const size_t k) {
    if(_binomialCache.size() > k)
        return _binomialCache[k];

    return _binomialCache.emplace_back(factorial(n) / (factorial(k) * factorial(n - k)));
}

BernsteinPolynomial::BernsteinPolynomial(const std::function<double(double)>& function, size_t degree) : _function {function}, _degree {degree} {
    _coefficients = computeCoefficients();
    _binomialCache = {};
}

std::vector<ExactDouble> BernsteinPolynomial::computeCoefficients() {
    std::vector<ExactDouble> coefficients = std::vector<ExactDouble>(_degree + 1);

    for(size_t i = 0; i <= _degree; ++i)
        coefficients[i] = ExactDouble(_function(i / (float)(_degree)));

    return coefficients;
}

Bounds<FloatMP> BernsteinPolynomial::bernsteinBasisPolynomialFor(size_t v, double x) {
    return binomialCoefficient(_degree, v) * ExactDouble(pow(x, v)) * ExactDouble(pow(1 - x, _degree - v));
}

Bounds<FloatMP> BernsteinPolynomial::evaluate(double x) {
    if(x < 0 || x > 1) return Bounds<FloatMP>(FloatMP(-inf, MultiplePrecision(128)), FloatMP(inf, MultiplePrecision(128)));

    Ariadne::Bounds<FloatMP> sum = FloatMP(0, MultiplePrecision(128));

    for(size_t i = 0; i <= _degree; ++i) {
        auto bp = bernsteinBasisPolynomialFor(i, x);
        sum += _coefficients[i] * bp;
    }

    return sum;
}

std::vector<Bounds<FloatMP>> BernsteinPolynomial::_factorialCache = {FloatMP(1, MultiplePrecision(128)), FloatMP(1, MultiplePrecision(128))};

}
