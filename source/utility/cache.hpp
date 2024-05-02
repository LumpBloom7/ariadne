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

// An alternative version that allows for initialization using passed function without the need of creating a subclass of Cache
template<typename ReturnT, typename... ArgT>
class SimpleCache : public Cache<ReturnT, ArgT...> {
  public:
    SimpleCache(std::function<ReturnT(ArgT...)> function) : _function{function} {}

    template<typename FT>
    SimpleCache(FT function) : _function{function} {}

  protected:
    ReturnT compute(ArgT... args) { return _function(args...); }

  private:
    std::function<ReturnT(ArgT...)> _function;
};
} // namespace Ariadne

#endif