#pragma once

#include "../Scheme.h"

namespace search_schemes::generator {

auto pigeon_trivial(int minK, int K) -> Scheme;
auto pigeon_opt(int minK, int K) -> Scheme;

}
