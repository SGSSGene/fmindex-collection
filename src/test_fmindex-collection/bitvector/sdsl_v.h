// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "concepts.h"

#include <sdsl/bit_vectors.hpp>

struct SDSL_V {
    sdsl::bit_vector bitvector;
    sdsl::bit_vector::rank_1_type bv;

    SDSL_V() = default;
    SDSL_V(SDSL_V&&) = default;
    SDSL_V(std::vector<bool> input)
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
    void serialize(Archive& ar) {
        ar(bitvector, bv);
    }
};
