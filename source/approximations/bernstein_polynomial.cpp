#include "approximations/bernstein_polynomial.hpp"

#include "numeric/rational.hpp"
#include "utility/hash_numeric.hpp"
#include "utility/hash_tuple.hpp"

namespace Ariadne {
SimpleCache<Integer, Nat64, Nat64> IBernsteinPolynomialBase::binomialCoefficients{[](auto n, auto m) {
    auto a = Factorials::get(n);
    auto b = Factorials::get(m);
    auto c = Factorials::get(n - m);

    return (a / (b * c)).get_num(); // This will always be an integer
}};
} // namespace Ariadne
