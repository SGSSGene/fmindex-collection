#pragma once

#include "../Scheme.h"

namespace oss::generator {

auto pigeon(int minK, int K, int sigma) -> Scheme;
auto pigeon_trivial(int minK, int K) -> Scheme;
auto pigeon_opt(int minK, int K) -> Scheme;

}
