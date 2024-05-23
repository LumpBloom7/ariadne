#ifndef ARIADNE_BERNSTEIN_POLYNOMIAL_GENERATOR_HPP
#define ARIADNE_BERNSTEIN_POLYNOMIAL_GENERATOR_HPP

#include <exception>

#include "approximations/bounded_bernstein_polynomial.hpp"
#include "approximations/composite_bernstein_polynomial.hpp"

namespace Ariadne {

template<typename T>
auto createBernsteinPolynomialWith(const std::function<Bounds<T>(Bounds<T>)> &function, const PositiveUpperBound<T> epsilon) {
    DegreeType degree = 1;
    auto b = BoundedBernsteinPolynomial(function, epsilon, degree, epsilon.precision());

    while (!decide(b.maximumError() < epsilon)) {
        std::cout << degree << " " << b.maximumError() << std::endl;

        degree *= 2;
        b = BoundedBernsteinPolynomial(function, epsilon, degree, epsilon.precision());
    }

    return b;
}

template<typename T>
BoundedBernsteinPolynomial<T> createIterativeBernsteinPolynomialWith(const std::function<Bounds<T>(Bounds<T>)> &function, const PositiveUpperBound<T> epsilon, int maxIterations = 10) {
    DegreeType degree = 1;

    while (true) {
        std::cout << "Degree: " << degree << std::endl;
        auto composite = BoundedBernsteinPolynomial(function, epsilon, degree, epsilon.precision());
        for (int i = 1; i < maxIterations; ++i) {
            std::cout << "\tMaximum error: " << composite.maximumError() << std::endl;
            std::cout << "\tIteration: " << i << std::endl;
            if (decide(composite.maximumError() < epsilon))
                return composite;

            std::function<Bounds<T>(Bounds<T>)> f = [&](Bounds<T> x) {
                return function(x) - composite(x);
            };

            auto p = BernsteinPolynomial<T>(f, degree, epsilon.precision());
            composite += p;

            // Recompute composite
            composite.computeErrorBounds(function, epsilon);
        }

        degree *= 2;
    }

    throw std::runtime_error("Coudn't compute iterative Bernstein via iteration");
}

} // namespace Ariadne

#endif
