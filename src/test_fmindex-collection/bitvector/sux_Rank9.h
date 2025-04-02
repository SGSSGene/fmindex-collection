// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "concepts.h"

#include <sdsl/bit_vectors.hpp>
#include <rank9.h>

struct Rank9 {
    std::vector<uint64_t> bitvector;
    rank9 bv;

    Rank9() = default;
    Rank9(Rank9&&) = default;
    Rank9(std::vector<bool> input)
        : bitvector{[&]() {
            auto bitvector = std::vector<uint64_t>{};
            bitvector.resize(input.size()/64 + 1);
            for (size_t i = 0; i < input.size(); ++i) {
                auto id = i / 64;
                auto offset = i % 64;
                bitvector[id] |= (input[i]<<offset);
            }
            return bitvector;
        }()}
        , bv{bitvector.data(), input.size()+1}
    {}

    auto symbol(size_t i) const -> bool {
        auto id = i / 64;
        auto offset = i % 64;
        return (bitvector[id] >> offset) & 1;
    }

    auto rank(size_t i) const -> size_t {
        //hack, since rank9 is not const correct
        return const_cast<rank9&>(bv).rank(i);
    }
    auto space_usage() const -> size_t {
        return 0;
//        bitvector.size() * 8 + rank9.
    }
};
