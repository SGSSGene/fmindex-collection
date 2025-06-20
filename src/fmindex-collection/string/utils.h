// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "concepts.h"

namespace fmc::string {

template <String_c String>
auto computeAccumulatedC(String const& str) -> std::array<size_t, String::Sigma+1> {
    auto C = std::array<size_t, String::Sigma+1>{0};
    for (size_t i{0}; i < String::Sigma; ++i) {
        C[i+1] = str.rank(str.size(), i) + C[i];
    }
    return C;
}

}
