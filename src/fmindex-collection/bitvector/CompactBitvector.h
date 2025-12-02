// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "concepts.h"
#include "../bitset_popcount.h"

#include <bit>
#include <bitset>
#include <cassert>
#include <ranges>
#include <vector>

namespace fmc::bitvector {

/**
 * CompactBitvector with interleaved superblocks, blocks and bits
 *
 * - Each group consist of 384bits, divided into 6 blocks.
 * - Each block uses 9bits to represents a value (6*9bits = 54bits).
 *   The last 10bits are padding bits, not used for any thing.
 * - Superblock consist of a single 64bit number
 *
 *   For 384bits, we need 512bits, or 1.333bits to save a single bit
 */
struct CompactBitvector {
    struct alignas(64) Superblock {
        uint64_t superBlockEntry{};
        uint64_t blockEntries{};
        std::array<uint64_t, 6> bits{};

        uint64_t rank(size_t idx) const noexcept {
            assert(idx < 384);

            auto blockId = idx >> 6;
            auto block = getBlock(blockId);
            auto keep = (uint64_t{idx} & 63);
            if (keep == 0) return superBlockEntry + block;

            auto maskedBits = bits[blockId] << (64 - keep);
            auto bitcount   = std::bitset<64>{maskedBits}.count();

            auto total = superBlockEntry + block + bitcount;
            return total;
        }

        bool symbol(size_t idx) const noexcept {
            assert(idx < 384);

            auto blockId = idx >> 6;
            auto bitId = idx & 63;
            return bits[blockId] & (uint64_t{1} << bitId);
        }

        void setBlock(size_t blockId, size_t value) {
            blockEntries = blockEntries & ~(uint64_t{0b111111111} << blockId*9);
            blockEntries = blockEntries | (uint64_t{value} << blockId*9);
        }
        auto getBlock(size_t blockId) const -> size_t{
            return uint64_t{0b111111111} & (blockEntries >> (blockId * 9));
        }


        struct Proxy {
            Superblock& sb;
            size_t      bid;
            auto operator=(size_t v) {
                sb.setBlock(bid, v);
            }

            friend auto operator+(Proxy const& p, size_t v) -> size_t {
                return v + p.sb.getBlock(p.bid);
            }
            friend auto operator+(size_t v, Proxy const& p) -> size_t {
                return v + p.sb.getBlock(p.bid);
            }

        };
        auto block(size_t blockId) {
            return Proxy{*this, blockId};
        }

        template <typename Archive>
        void serialize(this auto&& self, Archive& ar) {
            ar(self.superBlockEntry, self.blockEntries, self.bits);
        }
    };

    static constexpr size_t Sigma = 2;

    std::vector<Superblock> superblocks{Superblock{}};
    size_t                  totalLength{};

    template <typename CB>
    struct Layer {
        CB get;

        auto operator[](size_t i) -> decltype(auto) {
            return get(i);
        }
    };

    template <typename CB>
    CompactBitvector(size_t length, CB cb)
        : CompactBitvector{std::views::iota(size_t{}, length) | std::views::transform([&](size_t i) {
            return cb(i);
        })}
    {}

    template <std::ranges::sized_range range_t>
        requires std::convertible_to<std::ranges::range_value_t<range_t>, uint8_t>
    CompactBitvector(range_t&& _range) {

        reserve(_range.size());

        auto _length = _range.size();
        superblocks.resize(_length/(64*6) + 1);

        auto l0 = Layer{[&](size_t i) -> uint64_t& {
            return superblocks[i].superBlockEntry;
        }};
        auto l1 = Layer{[&](size_t i) {
            auto sbid = i / 6;
            auto& sb = superblocks[sbid];

            auto bid = i % 6;
            return sb.block(bid);
        }};
        auto l2 = Layer{[&](size_t i) -> uint64_t& {
            auto sbid = i / 6;
            auto& sb = superblocks[sbid];

            auto bid = i % 6;
            return sb.bits[bid];
        }};

//        auto acc = Accumulator<64, 384>{};

        for (auto iter = _range.begin(); iter != _range.end(); ++iter) {
            // run only if full block
            auto restBits = std::min(size_t{64}, _range.size() - totalLength);

            // concatenate next full block
            uint64_t bits = *iter;
            for (size_t i{1}; i < restBits; ++i) {
                bool value = *(++iter);
                bits = bits | (uint64_t{value} << i);
            }

            // update bits and blocks
            auto l0_id = totalLength / 384;
            auto l1_id = totalLength / 64;
            l2[l1_id] = bits;
            totalLength += restBits;

            // abort - if block not full,
            if (restBits < 64) {
                break;
            }
            l1[l1_id+1] = l1[l1_id] + std::popcount(bits);

            // check if next superblock is full
            if (totalLength % 384 == 0) {
                l0[l0_id+1] = l0[l0_id] + l1[l1_id+1];
                l1[l1_id+1] = 0;
            }
        }
    }

    CompactBitvector() = default;
    CompactBitvector(CompactBitvector const&) = default;
    CompactBitvector(CompactBitvector&&) noexcept = default;
    auto operator=(CompactBitvector const&) -> CompactBitvector& = default;
    auto operator=(CompactBitvector&&) noexcept -> CompactBitvector& = default;

    void reserve(size_t _length) {
        superblocks.reserve(_length/(64*6) + 1);
    }

    void push_back(bool _value) {
        auto blockId = (totalLength >>  6) % 6;
        auto bitId   = totalLength % 64;
        auto& bits   = superblocks.back().bits[blockId];
        bits         = bits | (size_t{_value} << bitId);
        totalLength += 1;

        if (totalLength % 64 == 0) { // new block
            auto newS = superblocks.back().getBlock(blockId) + std::popcount(superblocks.back().bits[blockId]);
            if (totalLength % 384 == 0) { // new super block + new block
                auto s = superblocks.back().superBlockEntry;
                superblocks.emplace_back();
                superblocks.back().superBlockEntry = s + newS;
            } else {
                superblocks.back().setBlock(blockId + 1, newS);
            }
        }
    }


    size_t size() const noexcept {
        return totalLength;
    }

    bool symbol(size_t idx) const noexcept {
        auto superblockId = idx / 384;
        auto bitId        = idx % 384;
        return superblocks[superblockId].symbol(bitId);
    }

    uint64_t rank(size_t idx) const noexcept {
        assert(idx <= size());
        auto superblockId = idx / 384;
        auto bitId        = idx % 384;
        auto v = superblocks[superblockId].rank(bitId);
        assert(v <= idx);
        return v;
    }



    template <typename Archive>
    void serialize(this auto&& self, Archive& ar) {
        ar(self.totalLength, self.superblocks);
    }
};
static_assert(Bitvector_c<CompactBitvector>);

}
