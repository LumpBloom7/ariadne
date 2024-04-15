#ifndef ARIADNE_HASH_NUMERIC
#define ARIADNE_HASH_NUMERIC

#include "numeric/integer.hpp"

namespace Ariadne {
template<>
struct Hash<Nat32> {
    size_t operator()(const Nat32& x) const { return Hash<uint32_t>()(x); }
};

template<>
struct Hash<Nat64> {
    size_t operator()(const Nat64& x) const { return Hash<uint64_t>()(x); }
};

template<>
struct Hash<Int32> {
    size_t operator()(const Int32& x) const { return Hash<int32_t>()(x); }
};

template<>
struct Hash<Int64> {
    size_t operator()(const Int64& x) const { return Hash<int64_t>()(x); }
};
} // namespace Ariadne

#endif