// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../bitset_popcount.h"
#include "concepts.h"
#include "utils.h"
#include "EPRV3.h"

#include <bit>
#include <vector>

#if __has_include(<cereal/types/bitset.hpp>)
#include <cereal/types/bitset.hpp>
#endif

#include <iostream>

#include "../ternarylogic.h"

#include "NEPRV8.h"

namespace fmindex_collection::string {


template <size_t TSigma, size_t l1_bits_ct, size_t l0_bits_ct, bool Align=true>
struct L0L1_NEPRV9 {
    static_assert(l1_bits_ct < l0_bits_ct, "first level must be smaller than second level");
    static_assert(l0_bits_ct-l1_bits_ct <= std::numeric_limits<uint16_t>::max(), "l0_bits_ct can only hold up to uint16_t bits");

    static constexpr size_t Sigma = TSigma;

    // number of full length bit vectors needed `2^bitct > TSigma`
    static constexpr auto bitct = std::bit_width(TSigma-1);
    // next full power of 2
    static constexpr auto bvct  = (1ull << bitct);

    struct alignas(Align?alignAsValue(bitct):size_t{1}) InBits {
        std::array<std::bitset<l1_bits_ct>, bitct> bits;

        uint8_t symbol(uint64_t idx) const {
            assert(idx < l1_bits_ct);
            uint8_t symb{};
            for (uint64_t i{bitct}; i > 0; --i) {
                auto b = bits[i-1].test(idx);
                symb = (symb<<1) | b;
            }
            return symb;
        }

        uint64_t rank(uint64_t idx, uint8_t symb) const {
            assert(idx <= l1_bits_ct);
            auto v = neprv8_detail::rank(bits, symb);
            return lshift_and_count(v, l1_bits_ct-idx);
        }

        uint64_t prefix_rank(uint64_t idx, uint64_t symb) const {
            assert(idx <= l1_bits_ct);
            auto v = neprv8_detail::prefix_rank(bits, symb);
            return lshift_and_count(v, l1_bits_ct-idx);
        }

        auto all_ranks(uint64_t idx) const -> std::array<uint64_t, TSigma> {
            assert(idx <= l1_bits_ct);

            auto vs = neprv8_detail::rank_all<(1ull<<bitct)>(bits);
            auto v = std::array<uint64_t, TSigma>{};
            static_assert(v.size() <= vs.size());
            for (size_t i{0}; i < v.size(); ++i) {
                auto count = skip_first_or_last_n_bits_and_count(vs[i], idx+l1_bits_ct);
                v[i] = count;
            }
            return v;
        }
        void setSymbol(size_t i, uint8_t symb) {
            assert(i <= l1_bits_ct);
            assert(symb < TSigma);
            for (size_t j{0}; j < bitct; ++j) {
                bits[j][i] = (symb>>j) & 1;
            }
        }

        template <typename Archive>
        void load(Archive& ar) {
            for (auto& v : bits) {
                loadBV(v, ar);
            }
        }
        template <typename Archive>
        void save(Archive& ar) const {
            for (auto const& v : bits) {
                saveBV(v, ar);
            }
        }
    };

    using BlockL1 = std::array<uint16_t, TSigma>;
    using BlockL0 = std::array<uint64_t, TSigma>;

    std::vector<InBits> bits{{}};
    std::vector<BlockL1> l1{{}};
    std::vector<BlockL0> l0{{}};
    size_t totalLength{};

    L0L1_NEPRV9() = default;
    L0L1_NEPRV9(std::span<uint8_t const> _symbols) {
        auto const _length = _symbols.size();
        bits.reserve(_length/l1_bits_ct + 2);
        if (_length == 0) return;

        // fill all inbits
        for (auto c : _symbols) {
            auto bitId         = totalLength % l1_bits_ct;
            bits.back().setSymbol(bitId, c);

            totalLength += 1;
            if (totalLength % l1_bits_ct == 0) { // next bit will require a new in-bits block
                bits.emplace_back();
            }
        }

        // fill l0/l1 structure
        {
            size_t l0BlockCt = (totalLength / l0_bits_ct) + 1;
            size_t l1BlockCt = l0BlockCt * (l0_bits_ct / l1_bits_ct);
            size_t inbitsCt  = l0BlockCt * ((l0_bits_ct) / l1_bits_ct);

            l0.resize(l0BlockCt);
            l1.resize(l1BlockCt);
            bits.resize(inbitsCt);

            constexpr size_t b1 = l1_bits_ct;
            constexpr size_t b0 = l0_bits_ct;

            constexpr size_t l1_block_ct = b0 / b1;


            BlockL0 l0_acc{};
            // walk through all superblocks
            for (size_t l0I{0}; l0I < l0BlockCt; ++l0I) {
                // update l0
                l0[l0I] = l0_acc;

                // right part
                BlockL1 acc{};
                for (size_t i{0}; i < l1_block_ct; ++i) {
                    auto idx = l0I*l1_block_ct + i;
                    l1[idx] = acc;
                    auto& b = bits[l0I*l1_block_ct + i];

                    auto counts = b.all_ranks(l1_bits_ct);

                    for (size_t symb{0}; symb < TSigma; ++symb) {
                        acc[symb] += counts[symb];
                    }
                }

                for (size_t symb{0}; symb < TSigma; ++symb) {
                    l0_acc[symb] += acc[symb];
                }
            }
        }
    }

    size_t size() const {
        return totalLength;
    }

    uint8_t symbol(uint64_t idx) const {
        assert(idx < totalLength);
        auto bitId = idx % l1_bits_ct;
        auto l1Id  = idx / l1_bits_ct;
        assert(l1Id < bits.size());

        return bits[l1Id].symbol(bitId);
    }

    uint64_t rank(uint64_t idx, uint8_t symb) const {
        assert(idx <= totalLength);
        auto bitId = idx % (l1_bits_ct);
        auto l1Id = idx / l1_bits_ct;
        auto l0Id = idx / l0_bits_ct;
        assert(l1Id < bits.size());
        assert(l0Id < l0.size());

        auto count = bits[l1Id].rank(bitId, symb);

        auto r =  l0[l0Id][symb] + l1[l1Id][symb] + count;
        assert(r <= idx);
        return r;
    }

    uint64_t prefix_rank(uint64_t idx, uint8_t symb) const {
        assert(idx <= totalLength);
        auto bitId = idx % (l1_bits_ct);
        auto l1Id = idx / l1_bits_ct;
        auto l0Id = idx / l0_bits_ct;
        assert(l1Id < bits.size());
        assert(l0Id < l0.size());

        size_t r{};
        for (size_t s{0}; s <= symb; ++s) {
            auto count = bits[l1Id].rank(bitId, s);
            r += l0[l0Id][s] + l1[l1Id][s] + count;
            assert(r <= idx);
        }
        return r;
    }

    auto all_ranks(uint64_t idx) const -> std::array<uint64_t, TSigma> {
        auto r = std::array<uint64_t, TSigma>{};
        for (size_t symb{0}; symb < TSigma; ++symb) {
            r[symb] = rank(idx, symb);
        }
        return r;
    }

    auto all_ranks_and_prefix_ranks(uint64_t idx) const -> std::tuple<std::array<uint64_t, TSigma>, std::array<uint64_t, TSigma>> {
        auto rs = all_ranks(idx);
        auto prs = rs;
        for (size_t i{1}; i < prs.size(); ++i) {
            prs[i] = prs[i] + prs[i-1];
        }
        return {rs, prs};
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(l0, l1, bits, totalLength);
    }
};

template <size_t Sigma> using L0L1_NEPRV9_64_4k   = L0L1_NEPRV9<Sigma, 64, 4096>;
template <size_t Sigma> using L0L1_NEPRV9_128_4k  = L0L1_NEPRV9<Sigma, 128, 4096>;
template <size_t Sigma> using L0L1_NEPRV9_256_4k  = L0L1_NEPRV9<Sigma, 256, 4096>;
template <size_t Sigma> using L0L1_NEPRV9_512_4k  = L0L1_NEPRV9<Sigma, 512, 4096>;
template <size_t Sigma> using L0L1_NEPRV9_1024_4k = L0L1_NEPRV9<Sigma, 1024, 4096>;
template <size_t Sigma> using L0L1_NEPRV9_2048_4k = L0L1_NEPRV9<Sigma, 2048, 4096>;

static_assert(checkRankVector<L0L1_NEPRV9_64_4k>);
static_assert(checkRankVector<L0L1_NEPRV9_128_4k>);
static_assert(checkRankVector<L0L1_NEPRV9_256_4k>);
static_assert(checkRankVector<L0L1_NEPRV9_512_4k>);
static_assert(checkRankVector<L0L1_NEPRV9_1024_4k>);
static_assert(checkRankVector<L0L1_NEPRV9_2048_4k>);

template <size_t Sigma> using L0L1_NEPRV9_64_64k   = L0L1_NEPRV9<Sigma, 64, 65536>;
template <size_t Sigma> using L0L1_NEPRV9_128_64k  = L0L1_NEPRV9<Sigma, 128, 65536>;
template <size_t Sigma> using L0L1_NEPRV9_256_64k  = L0L1_NEPRV9<Sigma, 256, 65536>;
template <size_t Sigma> using L0L1_NEPRV9_512_64k  = L0L1_NEPRV9<Sigma, 512, 65536>;
template <size_t Sigma> using L0L1_NEPRV9_1024_64k = L0L1_NEPRV9<Sigma, 1024, 65536>;
template <size_t Sigma> using L0L1_NEPRV9_2048_64k = L0L1_NEPRV9<Sigma, 2048, 65536>;

static_assert(checkRankVector<L0L1_NEPRV9_64_64k>);
static_assert(checkRankVector<L0L1_NEPRV9_128_64k>);
static_assert(checkRankVector<L0L1_NEPRV9_256_64k>);
static_assert(checkRankVector<L0L1_NEPRV9_512_64k>);
static_assert(checkRankVector<L0L1_NEPRV9_1024_64k>);
static_assert(checkRankVector<L0L1_NEPRV9_2048_64k>);

template <size_t Sigma> using L0L1_NEPRV9_64_64kUA   = L0L1_NEPRV9<Sigma, 64, 65536, false>;
template <size_t Sigma> using L0L1_NEPRV9_128_64kUA  = L0L1_NEPRV9<Sigma, 128, 65536, false>;
template <size_t Sigma> using L0L1_NEPRV9_256_64kUA  = L0L1_NEPRV9<Sigma, 256, 65536, false>;
template <size_t Sigma> using L0L1_NEPRV9_512_64kUA  = L0L1_NEPRV9<Sigma, 512, 65536, false>;
template <size_t Sigma> using L0L1_NEPRV9_1024_64kUA = L0L1_NEPRV9<Sigma, 1024, 65536, false>;
template <size_t Sigma> using L0L1_NEPRV9_2048_64kUA = L0L1_NEPRV9<Sigma, 2048, 65536, false>;

static_assert(checkRankVector<L0L1_NEPRV9_64_64kUA>);
static_assert(checkRankVector<L0L1_NEPRV9_128_64kUA>);
static_assert(checkRankVector<L0L1_NEPRV9_256_64kUA>);
static_assert(checkRankVector<L0L1_NEPRV9_512_64kUA>);
static_assert(checkRankVector<L0L1_NEPRV9_1024_64kUA>);
static_assert(checkRankVector<L0L1_NEPRV9_2048_64kUA>);

}
