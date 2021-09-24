#pragma once

#include "../Scheme.h"

namespace search_schemes::generator {

auto zeroOnesZero_trivial(int minK, int K) -> Scheme;
auto zeroOnesZero_opt(int minK, int K) -> Scheme;

}
