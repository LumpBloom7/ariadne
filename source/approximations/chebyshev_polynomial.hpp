#ifndef ARIADNE_CHEBYSHEV_POLYNOMIAL_HPP
#define ARIADNE_CHEBYSHEV_POLYNOMIAL_HPP
#include <cmath>
#include <functional>
#include <iostream>
#include <optional>

#include "mpfr.h"
#include "numeric/float_approximation.hpp"
#include "numeric/float_bounds.hpp"
#include "numeric/float_error.hpp"
#include "numeric/floatmp.hpp"
#include "numeric/integer.hpp"
#include "numeric/real.hpp"
#include "numeric/validated_real.hpp"
#include "utility/cache.hpp"
#include "utility/factorials.hpp"
#include "utility/hash.hpp"
#include "utility/hash_numeric.hpp"
#include "utility/hash_tuple.hpp"
#include "utility/standard.hpp"

namespace Ariadne {
template<typename T>
class ChebyshevPolynomial {
    typedef typename T::RoundingModeType RND;
    typedef typename T::PrecisionType PR;

  public:
    ChebyshevPolynomial(std::vector<Bounds<T>> coefficients) : _coefficients{coefficients} {}
    ChebyshevPolynomial(const std::function<Bounds<T>(Bounds<T>)>& function, int degree, size_t maxPrecision) {
        generateCoefficients(function, degree, MultiplePrecision(maxPrecision));
    }

    Bounds<T> evaluate(Bounds<T> x) {
        auto sum = -hlf(_coefficients[0]);

        for (int i = 0; i < _coefficients.size(); ++i)
            sum += _coefficients[i] * chebyshevPolynomial(x, i);

        return sum;
    }

  private:
    void generateCoefficients(const std::function<Bounds<T>(Bounds<T>)>& function, int degree, PR precisionType) {
        _coefficients = {};

        auto denominator = T(Approximation<T>(1 / static_cast<float>(degree), precisionType));
        auto factor = (2 * denominator);

        auto pi = T::pi(MPFR_RNDNA, precisionType);
        auto half = hlf(T(1, precisionType));

        for (int j = 1; j <= degree; ++j) {
            auto sum = Bounds<T>(precisionType);

            for (int k = 1; k <= degree; ++k) {
                auto kMinusHalf = k - half;

                auto a = cos(pi * kMinusHalf * denominator);
                auto b = cos(pi * (j - 1) * kMinusHalf * denominator);
                sum += function(a) * b;
            }

            _coefficients.emplace_back(factor * sum);
        }
    }

    // This is an alternate implementation of finding Tn that is much faster
    static Bounds<T> chebyshevPolynomialAlt(Bounds<T> x, int degree) {
        return cos(degree * acos(x));
    }

    static Bounds<T> chebyshevPolynomial(Bounds<T> x, int degree) {
        auto T0 = Bounds<T>(1, x.precision());

        if (degree == 0)
            return T0;

        auto T1 = x;

        if (degree == 1)
            return T1;

        auto TwoX = 2 * x;

        for (int i = 1; i < degree; ++i) {
            auto T2 = (TwoX * T1) - T0;

            T0 = T1;
            T1 = T2;
        }

        return T1;
    }

    std::vector<Bounds<T>> _coefficients;
};

} // namespace Ariadne

#endif
