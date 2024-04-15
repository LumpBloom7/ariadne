#ifndef ARIADNE_CACHE
#define ARIADNE_CACHE

#include <tuple>
#include <unordered_map>
#include "hash.hpp"

namespace Ariadne {
template<typename ReturnT, typename... ArgT>
class Cache {
  public:
    ReturnT operator()(ArgT... args) { return get(args...); };
    ReturnT get(ArgT... args) {
        auto retIt = _backingStorage.find(std::make_tuple(args...));
        if (retIt != _backingStorage.end()) {
            return retIt->second;
        }

        auto computed = compute(args...);

        _backingStorage.emplace(std::make_tuple(args...), computed);

        return computed;
    }

  protected:
    virtual ReturnT compute(ArgT... args) = 0;

  private:
    std::unordered_map<std::tuple<ArgT...>, ReturnT, Hash<std::tuple<ArgT...>>> _backingStorage = {};
};
} // namespace Ariadne

#endif
