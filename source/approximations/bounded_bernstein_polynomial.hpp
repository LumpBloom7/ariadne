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
    BoundedBernsteinPolynomial(const std::function<Bounds<T>(Bounds<T>)> &function, DegreeType degree, PR precision)
        : BernsteinPolynomial<T>(function, degree, precision) {
        test2(function, PositiveUpperBound<T>(T::inf(precision)));
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
            degree *= 2;

        } while (/* !endpointTest(function, targetEpsilon) || */ !test2(function, targetEpsilon));
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

    bool test2(const std::function<Bounds<T>(Bounds<T>)> &function, const PositiveUpperBound<T> &targetEpsilon) {
        auto degree = (this->_coefficients).size() - 1;

        _errorBounds.clear();
        _errorBounds.reserve(degree);

        auto denominator = invert(T(degree, targetEpsilon.precision()));

        auto rightXBounds = T(1, targetEpsilon.precision());
        Bounds<T> leftX = T(targetEpsilon.precision());
        Bounds<T> rightX = leftX + denominator;
        auto leftVal = this->evaluate(leftX);
        auto rightVal = this->evaluate(rightX);

        for (int i = 1; i <= degree; ++i) {
            auto minimum = min(leftVal.lower_raw(), rightVal.lower_raw());
            auto maximum = max(leftVal.upper_raw(), rightVal.upper_raw());

            auto domain = Bounds<T>(LowerBound<T>(leftX), UpperBound<T>(min(rightX, rightXBounds)));
            auto actual = function(domain);
            auto predicted = Bounds<T>(LowerBound<T>(minimum), UpperBound<T>(maximum));

            auto errorBounds = actual - predicted;

            auto errorUpperBound = mag(errorBounds);

            if ((errorUpperBound > targetEpsilon).repr() >= LogicalValue::INDETERMINATE)
                return false;

            _errorBounds.emplace_back(errorUpperBound);

            // Prep for next iteration
            leftX = rightX;
            rightX += denominator;

            leftVal = rightVal;
            rightVal = this->evaluate(rightX);
        }

        return true;
    }

    bool test(const std::function<Bounds<T>(Bounds<T>)> &function, const PositiveUpperBound<T> &targetEpsilon) {
        auto degree = (this->_coefficients).size() - 1;

        std::deque<PositiveUpperBound<T>> queue = {};
        auto denominator = invert(T(degree, targetEpsilon.precision()));

        auto intervalDomain = Bounds<T>(LowerBound<T>(targetEpsilon.precision()), UpperBound<T>(denominator));
        auto centre = (degree / 2 + 1);
        auto centreX = centre * denominator;

        // I evaluate error bounds from the center of our domain first
        // This is to take advantage of the fact that the BernsteinPolynomial tends to

        for (int i = 1; i < degree / 2 + 1; ++i) {
            for (int j = 0; j < 2; ++j) {
                int factor = 1 - j * 2;
                auto iF = (i * factor);

                if ((iF + centre) > degree || (iF + centre) < 0)
                    continue;

                auto x = centreX + (iF * denominator + intervalDomain);

                auto actual = function(x);
                auto predicted = this->evaluate(x);

                auto maxError = mag((actual - predicted) * 5 / 4);
                /*
                                std::cout << x << std::endl;
                                std::cout << iF << std::endl;
                                std::cout << actual << std::endl;
                                std::cout << predicted << std::endl;
                                std::cout << maxError << std::endl
                                          << degree << std::endl
                                          << std::endl;
                 */
                if ((maxError > targetEpsilon).repr() >= LogicalValue::INDETERMINATE)
                    return false;

                if (factor > 0)
                    queue.emplace_back(maxError);
                else
                    queue.emplace_front(maxError);
            }
        }

        _errorBounds.clear();
        _errorBounds.reserve(degree);
        while (!queue.empty()) {
            _errorBounds.emplace_back(queue.front());
            queue.pop_front();
        };

        return true;
    }

    bool computeErrorBounds(const std::function<Bounds<T>(Bounds<T>)> &function, const PositiveUpperBound<T> &targetEpsilon) {
        auto degree = (this->_coefficients).size() - 1;
        _errorBounds.clear();
        _errorBounds.reserve(degree);

        auto denominator = invert(T(degree, targetEpsilon.precision()));
        auto fDomain = Bounds<T>(LowerBound<T>(0, targetEpsilon.precision()), UpperBound<T>(1, targetEpsilon.precision()));
        auto fRange = function(fDomain);

        for (int i = 1; i <= degree; ++i) {
            auto x = i * denominator;
            auto interval = Bounds<T>(LowerBound<T>((i - 1) * denominator), UpperBound<T>(x));

            auto originalBounds = function(interval);
            auto polynomialBounds = this->evaluate(interval);

            // polynomialBounds = refinement(polynomialBounds, fRange);

            auto maxError = mag((originalBounds - polynomialBounds) * 5 / 4);

            std::cout << interval << std::endl;
            std::cout << originalBounds << std::endl;
            std::cout << polynomialBounds << std::endl;
            std::cout << maxError << std::endl
                      << originalBounds - polynomialBounds << std::endl
                      << degree << std::endl
                      << std::endl;

            if ((maxError > targetEpsilon).repr() >= LogicalValue::LIKELY)
                return false;

            _errorBounds.emplace_back(maxError);
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
