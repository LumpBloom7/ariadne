#ifndef ARIADNE_HASH
#define ARIADNE_HASH

namespace Ariadne {

// Forward non-specialized templates to std::hash
template<typename T>
struct Hash {
    size_t operator()(const T& t) const {
        return std::hash<T>()(t);
    }
};
} // namespace Ariadne

#endif
