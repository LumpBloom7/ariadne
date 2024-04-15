#ifndef ARIADNE_HASH_GMP
#define ARIADNE_HASH_GMP

#include <functional>

#include <gmpxx.h>

#include "hash.hpp"

// http://daniel-strecker.com/blog/2020-07-06_c++_hash_for_gmp_mpz_t_mpz_class/
// Fill in the license later 

namespace Ariadne {
template<>
struct Hash<mpz_srcptr> {
    size_t operator()(const mpz_srcptr x) const;
};

template<>
struct Hash<mpz_t> {
    size_t operator()(const mpz_t x) const;
};

template<>
struct Hash<mpz_class> {
    size_t operator()(const mpz_class& x) const;
};
} // namespace Ariadne
#endif
