#ifndef ARIADNE_FACTORIALS
#define ARIADNE_FACTORIALS

#include <vector>

#include "numeric/integer.hpp"

namespace Ariadne {
class Factorials {
  public:
    static const Integer& get(Nat64 x);

  private:
    static const Integer& factorial(Nat64 x);

    static std::vector<Integer> _backingStorage;
};
} // namespace Ariadne

#endif
