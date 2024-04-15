#ifndef ARIADNE_FACTORIALS
#define ARIADNE_FACTORIALS

#include <vector>

#include "numeric/integer.hpp"

namespace Ariadne {
class Factorials {
  public:
    static Integer get(Nat64 x);

  private:
    static Integer factorial(Nat64 x);

    static std::vector<Integer> _backingStorage;
};
} // namespace Ariadne

#endif
