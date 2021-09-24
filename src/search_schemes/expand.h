#pragma once

#include "Scheme.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>
#include <optional>
#include <type_traits>

namespace search_schemes {

auto expandCount(size_t oldLen, size_t newLen) -> std::vector<size_t>;

bool isExpandable(Search s, size_t newLen);

auto expandPI(std::vector<size_t> const& pi, size_t _newLen) -> std::vector<size_t>;

auto expandLowerBound(std::vector<size_t> const& pi, std::vector<size_t> bound, size_t _newLen) -> std::vector<size_t>;

auto expandUpperBound(std::vector<size_t> const& pi, std::vector<size_t> bound, size_t _newLen) -> std::vector<size_t>;

auto expand(Search s, size_t newLen) -> std::optional<Search>;
auto expand(Scheme ss, size_t newLen) -> Scheme;

}
