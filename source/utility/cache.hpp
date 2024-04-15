#ifndef ARIADNE_CACHE
#define ARIADNE_CACHE

// http://daniel-strecker.com/blog/2020-07-06_c++_hash_for_gmp_mpz_t_mpz_class/
// Fill in the license later 

#include <unordered_map>
#include "hash.hpp"

namespace Ariadne {
template<typename ReturnT, typename ArgT>
class Cache {
  public:
    ReturnT get(ArgT args) {
        auto retIt = _backingStorage.find(args);
        if (retIt != _backingStorage.end()) {
            return retIt->second;
        }

        auto computed = compute(args);

        _backingStorage.emplace(args, computed);

        return computed;
    }

  protected:
    virtual ReturnT compute(ArgT args) = 0;

  private:
    std::unordered_map<ArgT, ReturnT, Hash<ArgT>> _backingStorage = std::unordered_map<ArgT, ReturnT, Hash<ArgT>>();
};
} // namespace Ariadne

#endif
