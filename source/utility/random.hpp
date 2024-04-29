#ifndef ARIADNE_RANDOM_HPP
#define ARIADNE_RANDOM_HPP

#include <cstddef>

#include <gmpxx.h>

#include "numeric/floatmp.hpp"

namespace Ariadne {

FloatMP random(MultiplePrecision precision) {
    FloatMP tmp = FloatMP(pr);
    mpfr_urandom(static_cast<mpfr_ptr>(tmp.get_mpfr()), s, MPFR_RNDNA);
    return tmp;
}

} // namespace Ariadne

#endif
