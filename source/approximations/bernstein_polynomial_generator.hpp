#ifndef ARIADNE_BERNSTEIN_POLYNOMIAL_GENERATOR_HPP
#define ARIADNE_BERNSTEIN_POLYNOMIAL_GENERATOR_HPP

#include <exception>

#include "approximations/bounded_bernstein_polynomial.hpp"
#include "approximations/composite_bernstein_polynomial.hpp"

namespace Ariadne {

template<typename T>
BoundedBernsteinPolynomial<T> createBernsteinPolynomialWith(const std::function<Bounds<T>(Bounds<T>)> &function, const PositiveUpperBound<T> epsilon) {
    DegreeType degree = 1;
    auto p = BernsteinPolynomial<T>(function, degree, epsilon.precision());
    auto b = BoundedBernsteinPolynomial<T>(p, function, epsilon);

    while (!decide(b.maximumError() < epsilon)) {
        std::cout << degree << " " << b.maximumError() << std::endl;

        degree *= 2;
        p = BernsteinPolynomial<T>(function, degree, epsilon.precision());
        b = BoundedBernsteinPolynomial<T>(p, function, epsilon);
    }

    return b;
}

template<typename T>
BoundedBernsteinPolynomial<T> createIterativeBernsteinPolynomialWith(const std::function<Bounds<T>(Bounds<T>)> &function, const PositiveUpperBound<T> epsilon, int maxIterations = 10) {
    DegreeType degree = 1;

    while (true) {
        std::cout << "Degree: " << degree << std::endl;
        CompositeBernsteinPolynomial<T> composite = {};
        for (int i = 0; i < maxIterations; ++i) {
            std::function<Bounds<T>(Bounds<T>)> f = [&](Bounds<T> x) {
                return function(x) - composite(x);
            };
            std::cout << "\tIteration: " << i << std::endl;
            auto p = BernsteinPolynomial<T>(f, degree, epsilon.precision());
            composite += p;

            auto b = BoundedBernsteinPolynomial<T>(composite, function, epsilon);
            std::cout << "\tMaximum error: " << b.maximumError() << std::endl;
            if (decide(b.maximumError() < epsilon))
                return b;
        }

        degree *= 2;
    }

    throw std::runtime_error("Coudn't compute iterative Bernstein via iteration");
}

} // namespace Ariadne

#endif
