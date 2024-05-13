#ifndef ARIADNE_BOUNDED_BERNSTEIN_POLYNOMIAL_HPP
#define ARIADNE_BOUNDED_BERNSTEIN_POLYNOMIAL_HPP

#include <cmath>
#include <deque>
#include <exception>
#include <functional>
#include <iostream>
#include <optional>
#include <stack>

#include "approximations/bernstein_polynomial.hpp"
namespace Ariadne {
template<typename T>
class BoundedBernsteinPolynomial : public BernsteinPolynomial<T> {
    typedef typename T::RoundingModeType RND;
    typedef typename T::PrecisionType PR;

  public:
    BoundedBernsteinPolynomial(const std::function<Bounds<T>(Bounds<T>)> &function, DegreeType degree, PR precision, int secantIterations = 5)
        : BernsteinPolynomial<T>(function, degree, precision, secantIterations) {
        computeErrorBounds(function, PositiveUpperBound<T>(T::inf(precision)));
    }

    BoundedBernsteinPolynomial(const std::function<Bounds<T>(Bounds<T>)> &function, const PositiveUpperBound<T> &targetEpsilon)
        : BernsteinPolynomial<T>({}) {
        DegreeType degree = 1;

        while (!earlyBoundsTest(function, degree, targetEpsilon)) {
            degree *= 2;
            std::cout << degree << std::endl;
        }

        std::clog << "EARLY PASS DONE" << std::endl;

        do {
            std::cout << degree << std::endl;
            if (degree == 0) {
                std::clog << "Unable to increase degree past 2^16, bailing out. Error bounds will be based on the previous degree tested." << std::endl;
                break;
            }

            this->generateCoefficients(function, degree, targetEpsilon.precision());
            this->findCriticalPoints(degree, targetEpsilon);
            degree *= 2;

        } while (/* !endpointTest(function, targetEpsilon) || */ !computeErrorBounds(function, targetEpsilon));
    }

    PositiveUpperBound<T> maximumError() const {
        auto maximum = _errorBounds[0];

        for (int i = 1; i < _errorBounds.size(); ++i) {
            maximum = max(maximum, _errorBounds[i]);
        }

        return maximum;
    }

    PositiveUpperBound<T> maximumErrorAt(const Bounds<T> &x) const {
        auto degree = (this->_coefficients).size() - 1;
        auto denominator = invert(T(degree, x.precision()));
        auto maximum = PositiveUpperBound<T>(x.precision());

        auto xr = 0 * denominator;

        for (size_t i = 0; i <= degree; ++i) {
            xr += denominator;

            if ((x.lower() >= xr).repr() >= LogicalValue::LIKELY)
                continue;

            maximum = max(maximum, _errorBounds[i]);

            if ((x.upper() <= xr).repr() >= LogicalValue::LIKELY)
                break;
        }

        return maximum;
    }

  private:
    bool earlyBoundsTest(const std::function<Bounds<T>(Bounds<T>)> &function, DegreeType degree, const PositiveUpperBound<T> &targetEpsilon) {
        auto denominator = invert(T(degree, targetEpsilon.precision()));
        for (int i = 1; i <= degree; ++i) {
            auto interval = Bounds<T>(LowerBound<T>((i - 1) * denominator), UpperBound<T>(i * denominator));

            auto range = function(interval);
            auto error = mag(range.upper().raw() - range.lower().raw());

            if ((error > targetEpsilon).repr() >= LogicalValue::INDETERMINATE)
                return false;
        }

        return true;
    }

    bool computeErrorBounds(const std::function<Bounds<T>(Bounds<T>)> &function, const PositiveUpperBound<T> &targetEpsilon) {
        auto degree = (this->_coefficients).size() - 1;

        _errorBounds.clear();
        _errorBounds.reserve(degree);

        auto denominator = invert(T(degree, targetEpsilon.precision()));

        auto x = Bounds<T>(LowerBound<T>(0, targetEpsilon.precision()), UpperBound<T>(denominator));

        for (int i = 1; i <= degree; ++i) {
            auto actual = function(x);
            auto predicted = this->evaluate(x);

            auto errorUpperBound = mag(actual-predicted);

            if ((errorUpperBound > targetEpsilon).repr() >= LogicalValue::INDETERMINATE)
                return false;

            _errorBounds.emplace_back(errorUpperBound);

            x += denominator;
        }

        return true;
    }

    static Bounds<T> invert(const Bounds<T> &x) {
        return 1 / x;
    }
    static Bounds<T> invert(const T &x) {
        return 1 / x;
    }

    std::vector<PositiveUpperBound<T>> _errorBounds{};
};

} // namespace Ariadne

#endif
