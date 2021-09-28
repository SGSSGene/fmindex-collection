#pragma once

#include "Scheme.h"

#include <algorithm>
#include <cassert>

namespace search_schemes {

/* checks if Scheme covers patterns up to a certain error
 *
 */
auto isComplete(Scheme const& ss, size_t minK, size_t maxK) -> bool;

}
