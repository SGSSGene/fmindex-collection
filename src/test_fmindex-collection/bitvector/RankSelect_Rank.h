// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "concepts.h"

#include <shared/rank.hpp>

template <size_t N=0>
struct RankSelect {
    using BV = std::variant_alternative_t<N, std::variant<
        shared::RankBasic32<shared::RankBasicType::RANK_BASIC_STANDARD>,
        shared::RankBasic32<shared::RankBasicType::RANK_BASIC_COMPRESSED_HEADERS>,
        shared::RankBasic64<shared::RankBasicType::RANK_BASIC_STANDARD>,
        shared::RankBasic64<shared::RankBasicType::RANK_BASIC_COMPRESSED_HEADERS>,
        shared::RankCF64,
        shared::RankMPE64<shared::RankMPEType::RANK_MPE1>,
        shared::RankMPE64<shared::RankMPEType::RANK_MPE2>,
        shared::RankMPE64<shared::RankMPEType::RANK_MPE3>
    >>;
    BV bv;
    size_t totalLength{};
    RankSelect() = default;
    RankSelect(RankSelect&&) = default;
    RankSelect(std::vector<bool> const& input)
        : totalLength{input.size()} {

        auto data = std::vector<uint8_t>{};
        data.resize((input.size()+7)/8);
        for (size_t i{}; i < input.size()/8; ++i) {
            uint8_t b =
                  (input[i*8+0] << 0)
                | (input[i*8+1] << 1)
                | (input[i*8+2] << 2)
                | (input[i*8+3] << 3)
                | (input[i*8+4] << 4)
                | (input[i*8+5] << 5)
                | (input[i*8+6] << 6)
                | (input[i*8+7] << 7);
            data[i] = b;
        }
        if (input.size() % 8 != 0) {
            uint8_t b{};
            for (size_t i{0}; i < input.size() % 8; ++i) {
                b = b | (data[input.size()/8 + i] << i);
            }
            data.back() = b;
        }
        bv.build(data.data(), data.size());
    }
    auto symbol(size_t i) const -> bool {
        assert(i < totalLength);
        //!TODO is there a more efficient way?
        return rank(i+1) - rank(i);
    }
    auto rank(size_t i) const -> size_t {
        assert(i <= totalLength);
        return const_cast<BV&>(bv).rank(i);
    }
    auto space_usage() const -> size_t {
        return const_cast<BV&>(bv).getSize();
    }
};
