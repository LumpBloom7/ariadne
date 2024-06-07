#ifndef ARIADNE_NEWTON_POLYNOMIAL_HPP
#define ARIADNE_NEWTON_POLYNOMIAL_HPP

#include "bernstein_polynomial_shared.hpp"
#include "numeric/float_approximation.hpp"
#include "numeric/float_bounds.hpp"
#include "numeric/float_error.hpp"
#include "numeric/floatmp.hpp"
#include "numeric/integer.hpp"
#include "numeric/real.hpp"
#include "numeric/validated_real.hpp"
#include "utility/hash.hpp"
#include "utility/hash_numeric.hpp"
#include "utility/hash_tuple.hpp"
#include "utility/standard.hpp"

namespace Ariadne {
template<typename T>
class NewtonPolynomial : virtual public IPolynomialApproximation<T> {
    using RND = T::RoundingModeType;
    using PR = T::PrecisionType;

  public:
    PR _precision;
    DegreeType _degree;

    NewtonPolynomial(std::vector<Bounds<T>> coefficients) : coefficients{coefficients},
                                                            _precision(coefficients[0].precision()), _degree(coefficients.size()), degreeReciprocal(rec(T(_degree, _precision))), zero(Bounds<T>(_precision)) {}
    NewtonPolynomial(const std::function<Bounds<T>(Bounds<T>)> &function, DegreeType degree, PR precision, int secantIters = 5)
        : _precision(precision), _degree(degree), degreeReciprocal(rec(T(degree, precision))), zero(Bounds<T>(precision)) {
        generateCoefficients(function);
    }

    virtual Bounds<T> evaluate(const Bounds<T> &x, int subInterval = 1) const override {
        Bounds<T> sum = zero;

        auto mult = zero + 1;

        for (int i = 0; i <= _degree; ++i) {
            sum += dividedDifferences[i][0] * mult;
            mult *= (x - xVals[i]);
        }

        return sum;
    }

    virtual Bounds<T> evaluateDerivative(const Bounds<T> &x) const override {
        Bounds<T> sum = zero;

        for (int i = 0; i < _degree; ++i) {
            Bounds<T> sum2 = zero;

            for (int j = 0; j < _degree; ++i) {
                if (j == i)
                    continue;

                sum2 += 1 / (x - xVals[j]);
            }

            auto mult1 = zero + 1;

            for (int j = 0; j < i; ++j) {
                mult1 *= (x - xVals[j]) / (xVals[i] - xVals[j]);
            }

            for (int j = i + 1; j < _degree; ++j) {
                mult1 *= (x - xVals[j]) / (xVals[i] - xVals[j]);
            }

            sum += mult1 * sum2 * coefficients[i];
        }
        return zero;
    }

    void generateCoefficients(const std::function<Bounds<T>(Bounds<T>)> &function) {
        coefficients.clear();
        coefficients.reserve(_degree + 1);

        auto x = zero;
        for (int i = 0; i <= _degree; ++i) {
            xVals.emplace_back(x);

            for (int j = 0; j <= i; ++j) {
                if (dividedDifferences.size() == j)
                    dividedDifferences.emplace_back(std::vector<Bounds<T>>{});

                if (j == 0)
                    dividedDifferences[j].emplace_back(function(xVals[i]));
                else {
                    auto top = dividedDifferences[j - 1][i - j];
                    auto bottom = dividedDifferences[j - 1][i - j + 1];

                    auto divisor = degreeReciprocal * j;
                    dividedDifferences[j].emplace_back((bottom - top) / divisor);
                }
            }

            x += degreeReciprocal;
        }
    }

    virtual DegreeType degree() const override {
        return _degree;
    }
    virtual PR precision() const override {
        return _precision;
    }

  private:
    std::vector<Bounds<T>> coefficients = {};
    std::vector<std::vector<Bounds<T>>> dividedDifferences = {};
    std::vector<Bounds<T>> xVals = {};
    Bounds<T> degreeReciprocal;
    Bounds<T> zero;
};
} // namespace Ariadne
#endif
