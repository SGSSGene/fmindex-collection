/*--

This file taken from libsais, a library for linear time suffix array,
longest common prefix array and burrows wheeler transform construction.

   Copyright (c) 2021-2022 Ilya Grebnov <ilya.grebnov@gmail.com>

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

Please see the file LICENSE for full copyright information.

--*/
#pragma once

#include "sais32.hpp"
#include "sais64.hpp"

#include <ranges>
#include <vector>

namespace libsais {

template <std::ranges::contiguous_range Range>
auto constructSA32(Range const range) -> std::vector<int32_t> {
    static_assert(std::is_same_v<uint8_t const*, decltype(range.data())>, "given range must be of type \"uint8_t const*\"");
    auto sa = std::vector<int32_t>(range.size());
    auto r = libsais32::constructSA(range.data(), sa.data(), range.size(), 0, nullptr);
    if (r != 0) { throw std::runtime_error("something went wrong constructing the SA"); }
    return sa;
}

template <std::ranges::contiguous_range Range>
auto constructSA64(Range const range) -> std::vector<int64_t> {
    static_assert(std::is_same_v<uint8_t const*, decltype(range.data())>, "given range must be of type \"uint8_t const*\"");
    auto sa = std::vector<int64_t>(range.size());
    auto r = libsais64::constructSA(range.data(), sa.data(), range.size(), 0, nullptr);
    if (r != 0) { throw std::runtime_error("something went wrong constructing the SA"); }
    return sa;
}

}
