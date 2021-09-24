#pragma once

#include "../Scheme.h"

namespace search_schemes::generator {

auto pigeon_trivial(size_t minK, size_t K) -> Scheme;
auto pigeon_opt(size_t minK, size_t K) -> Scheme;

}
