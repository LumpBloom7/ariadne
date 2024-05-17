#include "factorials.hpp"

namespace Ariadne {

const Integer& Factorials::get(Nat64 x) {
    if (x < _backingStorage.size())
        return _backingStorage[x];

    return factorial(x);
}

const Integer& Factorials::factorial(Nat64 x) {
    for (size_t i = _backingStorage.size(); i <= x; ++i) {
        _backingStorage.emplace_back(_backingStorage.back() * i);
    }

    return _backingStorage.back();
}

std::vector<Integer> Factorials::_backingStorage = {1, 1};

} // namespace Ariadne
