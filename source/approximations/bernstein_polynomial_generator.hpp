#ifndef ARIADNE_BERNSTEIN_POLYNOMIAL_GENERATOR_HPP
#define ARIADNE_BERNSTEIN_POLYNOMIAL_GENERATOR_HPP

#include <exception>
#include <memory>

#include "approximations/bounded_bernstein_polynomial.hpp"

namespace Ariadne {

template<typename T>
auto createBernsteinPolynomialWith(const std::function<Bounds<T>(Bounds<T>)> &function, const PositiveUpperBound<T> epsilon, int maxDepth = 1, int subIntervals = 1) {
    DegreeType degree = 1;

    std::shared_ptr<IPolynomialApproximation<T>> b = std::make_shared<BernsteinPolynomial<T>>(function, degree, epsilon.precision());
    auto bounded = BoundedPolynomialApproximation(function, b, epsilon, maxDepth, subIntervals);

    while (!decide(bounded.maximumError() < epsilon)) {
        std::cout << degree << " " << bounded.maximumError() << std::endl;

        degree *= 2;
        b = std::make_shared<BernsteinPolynomial<T>>(function, degree, epsilon.precision());
        bounded = BoundedPolynomialApproximation(function, b, epsilon, maxDepth, subIntervals);
    }

    return bounded;
}

template<typename T>
auto createIterativeBernsteinPolynomialWith(const std::function<Bounds<T>(Bounds<T>)> &function, const PositiveUpperBound<T> epsilon, int maxIterations = 10) {
    DegreeType degree = 1;

    while (true) {
        std::cout << "Degree: " << degree << std::endl;
        std::shared_ptr<IPolynomialApproximation<T>> composite = std::make_shared<BernsteinPolynomial<T>>(function, degree, epsilon.precision());
        auto boundedComposite = BoundedPolynomialApproximation(function, composite);

        for (int i = 1; i < maxIterations; ++i) {
            std::cout << "\tMaximum error: " << boundedComposite.maximumError() << std::endl;
            std::cout << "\tIteration: " << i << std::endl;
            if (decide(boundedComposite.maximumError() < epsilon))
                return boundedComposite;

            std::function<Bounds<T>(Bounds<T>)> f = [&](Bounds<T> x) {
                return function(x) - boundedComposite(x);
            };

            auto p = BernsteinPolynomial<T>(f, degree, epsilon.precision());
            *(dynamic_cast<BernsteinPolynomial<T> *>(composite.get())) += p;

            // Recompute composite
            boundedComposite.computeErrorBounds(function, epsilon);
        }

        degree *= 2;
    }

    throw std::runtime_error("Coudn't compute iterative Bernstein via iteration");
}

} // namespace Ariadne

#endif
