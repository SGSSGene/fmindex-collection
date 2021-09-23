#pragma once

#include "../Scheme.h"

namespace oss::generator {

auto zeroOnesZero(int minK, int K, int sigma, int DB_SIZE) -> Scheme;

auto zeroOnesZero_trivial(int minK, int K) -> Scheme;
auto zeroOnesZero_opt(int minK, int K) -> Scheme;


}
