#pragma once

#include "../Scheme.h"

namespace search_schemes::generator {

auto zeroOnesZero_trivial(size_t minK, size_t K) -> Scheme;
auto zeroOnesZero_opt(size_t minK, size_t K) -> Scheme;

}
