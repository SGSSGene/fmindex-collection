// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "concepts.h"

#include <sdsl/bit_vectors.hpp>

struct SDSL_V5 {
    sdsl::bit_vector bitvector;
    sdsl::rank_support_v5<> bv;

    SDSL_V5() = default;
    SDSL_V5(SDSL_V5&&) = default;
    SDSL_V5(std::vector<bool> input)
        : bitvector{[&]() {
            auto bv = sdsl::bit_vector(input.size(), 0);
            for (size_t i = 0; i < input.size(); ++i) {
                bv[i] = input[i];
            }
            return bv;
        }()}
        , bv{&bitvector}
    {}

    auto symbol(size_t i) const -> bool {
        return bitvector[i];
    }

    auto rank(size_t i) const -> size_t {
        return bv.rank(i);
    }

    template <typename Archive>
    void serialize(this auto&& self, Archive& ar) {
        ar(self.bitvector, self.bv);
    }
};
