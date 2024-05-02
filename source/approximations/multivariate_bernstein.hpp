#ifndef ARIADNE_MULTIVARIATE_BERNSTEIN
#define ARIADNE_MULTIVARIATE_BERNSTEIN

#include <array>
#include <cmath>
#include <cstddef>
#include <deque>
#include <exception>
#include <functional>
#include <iostream>
#include <optional>
#include <stack>
#include <tuple>

#include "bernstein_polynomial.hpp"

namespace Ariadne {
template<typename T, size_t argumentsNumber>
class MultivariateBernsteinPolynomial : BernsteinPolynomialBase {
    typedef typename T::RoundingModeType RND;
    typedef typename T::PrecisionType PR;

  public:
    MultivariateBernsteinPolynomial(const std::function<Bounds<T>(std::vector<Bounds<T>>)>& function, DegreeType degree, PR precision) : degree{degree}, denominator{T(1, precision) / T(degree, precision)} {
        generateCoefficients(function, degree, precision);
    }

    Bounds<T> evaluate(const std::vector<Bounds<T>>& args) const {
        std::array<int, argumentsNumber> indices = {};

        return evaluateRecursive(args, indices);
    }

    Bounds<T> operator()(const std::vector<Bounds<T>>& args) const { return evaluate(args); }

  private:
    Bounds<T> evaluateRecursive(const std::vector<Bounds<T>>& args, std::array<int, argumentsNumber>& indices, int level = 0) const {
        if (level == argumentsNumber) {
            auto actualIndex = flattenIndex(indices, degree + 1);
            return _coefficients[actualIndex];
        }

        Bounds<T> sum = Bounds<T>(args[level].precision());

        for (size_t i = 0; i <= degree; ++i) {
            indices[level] = i;
            auto bp = bernsteinBasisPolynomialFor(i, degree, args[level]);
            auto _coeff = evaluateRecursive(args, indices, level + 1);
            sum = fma(bp, _coeff, sum);
        }

        return sum;
    }

    void generateCoefficients(const std::function<Bounds<T>(std::vector<Bounds<T>>)>& function, DegreeType degree, PR precision) {
        _coefficients.clear();

        std::vector<Bounds<T>> arr = {};

        auto tmp = Bounds<T>(precision);
        for (int i = 0; i < argumentsNumber; ++i) {
            arr.emplace_back(tmp);
        }

        generateCoefficientsRecursive(function, degree, precision, arr);
    }

    void generateCoefficientsRecursive(const std::function<Bounds<T>(std::vector<Bounds<T>>)>& function, DegreeType degree, PR precision, std::vector<Bounds<T>>& array, int currentCount = 0) {
        if (currentCount == argumentsNumber) {
            auto res = function(array);
            _coefficients.emplace_back(res);
            return;
        }

        for (size_t i = 0; i <= degree; ++i) {
            auto x = i * denominator;

            array[currentCount] = x;
            generateCoefficientsRecursive(function, degree, precision, array, currentCount + 1);
        }
    }

    static Bounds<T> bernsteinBasisPolynomialFor(int v, int n, const Bounds<T>& x) {
        int v2 = (v * 2 >= n) ? (n - v) : v; // Used to avoid computing extra redundant values in binomial cache

        return binomialCoefficients(n, v2) * pow(x.value(), v) * pow(1 - x.value(), n - v);
    }

    static size_t flattenIndex(std::array<int, argumentsNumber> indices, size_t max) {
        size_t index = 0;
        for (int i = 0; i < argumentsNumber; ++i) {
            index += indices[i] * pow(max, (argumentsNumber - 1) - i);
        }
        return index;
    }

    std::vector<Bounds<T>> _coefficients = {};
    DegreeType degree;

    Bounds<T> denominator;
};
} // namespace Ariadne

#endif
