#include "factorials.hpp"

namespace Ariadne {

Integer Factorials::get(Nat64 x) {
    if (x < _backingStorage.size())
        return _backingStorage[x];

    return factorial(x);
}

Integer Factorials::factorial(Nat64 x) {
    auto res = _backingStorage.back();

    for (size_t i = _backingStorage.size(); i <= x; ++i) {
        res *= i;
        _backingStorage.emplace_back(res);
    }

    return res;
}

std::vector<Integer> Factorials::_backingStorage = {1, 1};

} // namespace Ariadne
