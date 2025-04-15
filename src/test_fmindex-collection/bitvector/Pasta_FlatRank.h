// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "concepts.h"

#include <pasta/bit_vector/bit_vector.hpp>
#include <pasta/bit_vector/support/rank.hpp>
#include <pasta/bit_vector/support/flat_rank.hpp>

struct FlatRank {
    pasta::BitVector bv;
    pasta::FlatRank<pasta::OptimizedFor::DONT_CARE, pasta::BitVector> rs{bv};
    FlatRank() = default;
    FlatRank(FlatRank&&) = default;
    FlatRank(std::vector<bool> const& input)
        : bv{[&]() {
            pasta::BitVector bv{input.size(), 0};
            for (size_t i = 0; i < bv.size(); ++i) {
                bv[i] = input[i];
            }
            return bv;
        }()}
    {
    }
    auto symbol(size_t i) const -> bool {
        return bv[i];
    }
    auto rank(size_t i) const -> size_t {
        return rs.rank1(i);
    }
    auto space_usage() const -> size_t {
        return bv.space_usage() + rs.space_usage();
    }
};
