#ifndef ARIADNE_CACHE
#define ARIADNE_CACHE

#include <tuple>
#include <unordered_map>

#include "hash.hpp"

namespace Ariadne {
template<typename ReturnT, typename... ArgT>
class Cache {
  public:
    const ReturnT& operator()(ArgT... args) { return get(args...); };
    const ReturnT& get(ArgT... args) {
        auto tpl = std::make_tuple(args...);
        auto retIt = _backingStorage.find(tpl);
        if (retIt != _backingStorage.end()) {
            return retIt->second;
        }

        return _backingStorage[tpl] = compute(args...);
    }

  protected:
    virtual const ReturnT compute(ArgT... args) = 0;

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
    const ReturnT compute(ArgT... args) { return _function(args...); }

  private:
    const std::function<ReturnT(ArgT...)> _function;
};
} // namespace Ariadne

#endif
