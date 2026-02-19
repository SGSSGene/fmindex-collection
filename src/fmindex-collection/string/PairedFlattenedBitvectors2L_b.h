// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../bitset_popcount.h"
#include "../ternarylogic.h"
#include "../utils.h"
#include "EPRV3.h"
#include "NEPRV8.h"
#include "concepts.h"
#include "utils.h"

#include <bit>
#include <limits>
#include <vector>

#if __has_include(<cereal/types/bitset.hpp>)
#include <cereal/types/bitset.hpp>
#endif

namespace fmc::string {

template <size_t TSigma, size_t L, size_t l1_bits_ct, size_t l0_bits_ct, bool Align=true>
struct PairedFlattenedBitvectors2L_b {
    static_assert(l1_bits_ct < l0_bits_ct, "first level must be smaller than second level");
    static_assert(l0_bits_ct-l1_bits_ct <= std::numeric_limits<uint16_t>::max(), "l0_bits_ct can only hold up to uint16_t bits");

    static constexpr size_t Sigma = TSigma;

    // number of full length bit vectors needed `2^bitct > TSigma`
    static constexpr auto bitct = std::bit_width(TSigma-1);
    // next full power of 2
    static constexpr auto bvct  = (1ull << bitct);

    static constexpr auto AlignVHigh = Align?alignAsValue(l1_bits_ct):alignof(std::array<std::bitset<l1_bits_ct>, L>);
    static constexpr auto AlignVLow  = Align?alignAsValue(l1_bits_ct):alignof(std::array<std::bitset<l1_bits_ct>, bitct<L?size_t{0}:(bitct-L)>);
    struct alignas(AlignVHigh) InBitsHigh {
        std::array<std::bitset<l1_bits_ct>, L> bits;

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
    struct alignas(AlignVLow) InBitsLow {
        std::array<std::bitset<l1_bits_ct>, bitct<L?size_t{0}:bitct-L> bits;

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

    struct InBits {
        std::array<std::bitset<l1_bits_ct>, bitct> bits;
    };


    using BlockL1 = std::array<uint16_t, TSigma>;
    using BlockL0 = std::array<uint64_t, TSigma>;

    mmser::vector<InBitsLow>  bits_low{{}};
    mmser::vector<InBitsHigh> bits_high{{}};
    mmser::vector<BlockL1> l1{{}};
    mmser::vector<BlockL0> l0{{}};


    size_t totalLength{};

    PairedFlattenedBitvectors2L_b()
        : PairedFlattenedBitvectors2L_b{internal_tag{}, std::span<uint8_t const>{}}
    {}

    PairedFlattenedBitvectors2L_b(std::span<uint8_t const> _symbols)
        : PairedFlattenedBitvectors2L_b{internal_tag{}, _symbols}
    {}

    PairedFlattenedBitvectors2L_b(std::span<uint64_t const> _symbols)
        : PairedFlattenedBitvectors2L_b{internal_tag{}, _symbols}
    {}

    template <std::ranges::range range_t>
        requires std::convertible_to<std::ranges::range_value_t<range_t>, uint64_t>
    PairedFlattenedBitvectors2L_b(range_t&& _symbols)
        : PairedFlattenedBitvectors2L_b{internal_tag{}, _symbols}
    {}


private:
    struct internal_tag{};

    template <std::ranges::range range_t>
        requires std::convertible_to<std::ranges::range_value_t<range_t>, uint64_t>
    PairedFlattenedBitvectors2L_b(internal_tag, range_t&& _symbols) {

        if constexpr (requires() { _symbols.size(); }) {
            auto const _length = _symbols.size();
            bits_low.reserve(_length/l1_bits_ct + 2);
            bits_high.reserve(_length/l1_bits_ct + 2);
        }

        // fill all inbits
        for (auto c : _symbols) {
            auto bitId = totalLength % l1_bits_ct;
            for (size_t j{0}; j < bitct; ++j) {
                if (j < bitct - L) {
                    bits_low.back().bits[j][bitId] = (c>>j) & 1;
                } else {
                    bits_high.back().bits[j - bitct + L][bitId] = (c>>j) & 1;
                }
            }

            totalLength += 1;
            if (totalLength % l1_bits_ct == 0) { // next bit will require a new in-bits block
                //bits.emplace_back();
                bits_low.emplace_back();
                bits_high.emplace_back();
            }
        }

        // fill l0/l1 structure
        {
            size_t l0BlockCt = (totalLength / (l0_bits_ct*2)) + 1;
            size_t l1BlockCt = l0BlockCt * (l0_bits_ct / l1_bits_ct);
            size_t inbitsCt  = l0BlockCt * ((l0_bits_ct*2) / l1_bits_ct);

            l0.resize(l0BlockCt);
            l1.resize(l1BlockCt);
            //bits.resize(inbitsCt);
            bits_low.resize(inbitsCt);
            bits_high.resize(inbitsCt);

            constexpr size_t b1 = l1_bits_ct;
            constexpr size_t b0 = l0_bits_ct;

            constexpr size_t l1_block_ct = b0 / b1;

            BlockL0 l0_acc{};
            // walk through all superblocks
            for (size_t l0I{0}; l0I < l0BlockCt; ++l0I) {
                // walk left to right and set l1 values (as if they are the beginning of a superblock)

                // left part
                BlockL0 acc{};
                for (size_t i{0}; i < l1_block_ct; ++i) {
                    auto idx = l0I*l1_block_ct*2 + i;
                    auto& b_low  = bits_low[idx];
                    auto& b_high = bits_high[idx];
                    auto counts = std::array<size_t, Sigma>{};
                    for (size_t symb{0}; symb < Sigma; ++symb) {
                        auto v = mark_rank_cb<l1_bits_ct, bitct>(symb, [&](size_t idx) -> auto const& {
                            if (idx < bitct - L) return b_low.bits[idx];
                            else                 return b_high.bits[idx-bitct+L];
                        });
                        counts[symb] = skip_first_or_last_n_bits_and_count(v, 0);
                    }
                    size_t a{};
                    for (size_t symb{0}; symb < TSigma; ++symb) {
                        a += counts[symb];
                        acc[symb] += a;
                    }
                    if (i % 2 == 0) {
                        for (size_t symb{0}; symb < TSigma; ++symb) {
                            l1[l0I*l1_block_ct + i/2][symb] = acc[symb];
                        }
                    }
                }
                for (size_t symb{0}; symb < TSigma; ++symb) {
                    l0_acc[symb] += acc[symb];
                }
                // update l0 (reached center)
                l0[l0I] = l0_acc;
                // walk backwards through left part and revert l0
                for (size_t i{0}; i < l1_block_ct; ++i) {
                    if (i % 2 == 0) {
                        for (size_t symb{0}; symb < TSigma; ++symb) {
                            auto idx = l0I*l1_block_ct + i/2;
                            l1[idx][symb] = acc[symb] - l1[idx][symb];
                        }
                    }
                }

                // right part
                acc = {};
                for (size_t i{l1_block_ct}; i < l1_block_ct*2; ++i) {
                    auto idx = l0I*l1_block_ct*2 + i;
                    auto& b_low  = bits_low[idx];
                    auto& b_high = bits_high[idx];
                    auto counts = std::array<size_t, Sigma>{};
                    for (size_t symb{0}; symb < Sigma; ++symb) {
                        auto v = mark_rank_cb<l1_bits_ct, bitct>(symb, [&](size_t idx) -> auto const& {
                            if (idx < bitct - L) return b_low.bits[idx];
                            else                 return b_high.bits[idx-bitct+L];
                        });
                        counts[symb] = skip_first_or_last_n_bits_and_count(v, 0);
                    }

                    size_t a{};
                    for (size_t symb{0}; symb < TSigma; ++symb) {
                        a += counts[symb];
                        acc[symb] += a;
                    }
                    if (i % 2 == 0) {
                        for (size_t symb{0}; symb < TSigma; ++symb) {
                            l1[l0I*l1_block_ct + i/2][symb] = acc[symb];
                        }
                    }
                }

                for (size_t symb{0}; symb < TSigma; ++symb) {
                    l0_acc[symb] += acc[symb];
                }
            }
        }
    }

public:
    size_t size() const {
        return totalLength;
    }

    uint64_t symbol(uint64_t idx) const {
        assert(idx < totalLength);
        auto bitId = idx % l1_bits_ct;
        auto l1Id  = idx / l1_bits_ct;
        assert(l1Id < bits_low.size());

        uint64_t symb{};
        for (uint64_t _i{bitct}; _i > 0; --_i) {
            auto i = _i-1;

            auto bit = [&]() {
                if (i < bitct - L) return bits_low[l1Id].bits[i].test(bitId);
                else               return bits_high[l1Id].bits[i-bitct+L].test(bitId);
            }();
            symb = (symb<<1) | bit;
        }
        assert(symb < Sigma);
        return symb;
    }

    template <size_t L2 = L>
    uint64_t symbol_limit(uint64_t idx) const {
        assert(idx < totalLength);
        auto bitId = idx % l1_bits_ct;
        auto l1Id  = idx / l1_bits_ct;
        assert(l1Id < bits_low.size());

        uint64_t symb{};
        for (uint64_t i{0}; i < L2; ++i) {
            auto bit = [&](size_t i) {
                i += (bitct - L2);
                if (i < bitct - L) return bits_low[l1Id].bits[i].test(bitId);
                else               return bits_high[l1Id].bits[i-bitct+L].test(bitId);
            }(i);
            symb = symb | (bit << i);
        }
        assert(symb < Sigma);
        return symb;
    }

    uint64_t rank(uint64_t idx, uint64_t symb) const {
        return rank_limit<bitct>(idx, symb);
    }

    template <size_t L2 = L>
    uint64_t rank_limit(uint64_t idx, uint64_t symb) const {
        assert(idx <= totalLength);
        assert(symb < (size_t{1}<<L2));
        auto bitId = idx % (l1_bits_ct*2);
        auto l1Id = idx / l1_bits_ct;
        auto l0Id = idx / l0_bits_ct;
        assert(l1Id < bits_low.size());
        assert(l1Id/2 < l1.size());
        assert(l0Id/2 < l0.size());

        auto count = [&]() {
            auto v = mark_rank_cb<l1_bits_ct, L2>(symb, [&](size_t idx) -> auto const& {
                idx += bitct - L2;
                if (idx < bitct - L) return bits_low[l1Id].bits[idx];
                else                 return bits_high[l1Id].bits[idx-bitct+L];
            });
            return skip_first_or_last_n_bits_and_count(v, bitId);
        }();

        auto symb1 = ((symb+1)<<(bitct-L2))-1;
        auto right_l1 = (l1Id%2);
        auto right_l0 = (l0Id%2);
        auto [superblock, block] = [&]() -> std::tuple<size_t, size_t> {
            if (symb == 0) return {l0[l0Id/2][symb1], l1[l1Id/2][symb1]};
            auto symb0 = (symb<<(bitct-L2))-1;
            return {
                (l0[l0Id/2][symb1] - l0[l0Id/2][symb0]),
                (l1[l1Id/2][symb1] - l1[l1Id/2][symb0])
            };
        }();
        auto r = [&]() {
            if (right_l0 && right_l1) {
                return superblock + block + count;
            } else if (right_l0) {
                return superblock + block - count;
            } else if (right_l1) {
                return superblock - block + count;
            } else {
                return superblock - block - count;
            }
        }();
        assert(r <= idx);
        return r;
    }


    uint64_t prefix_rank(uint64_t idx, uint64_t symb) const {
        return prefix_rank_limit<bitct>(idx, symb);
    }

    template <size_t L2 = L>
    uint64_t prefix_rank_limit(uint64_t idx, uint64_t symb) const {
        if (symb == 0) return 0;
        if (symb == (size_t{1}<<L2)) return idx;
        assert(idx <= totalLength);
        assert(symb <= (size_t{1}<<L2));
        auto bitId = idx % (l1_bits_ct*2);
        auto l1Id = idx / l1_bits_ct;
        auto l0Id = idx / l0_bits_ct;
        assert(l1Id < bits_low.size());
        assert(l1Id/2 < l1.size());
        assert(l0Id/2 < l0.size());

        auto count = [&]() {
            auto v = mark_prefix_rank_cb<l1_bits_ct, L2>(symb, [&](size_t idx) -> auto const& {
                idx += bitct - L2;
                if (idx < bitct - L) return bits_low[l1Id].bits[idx];
                else                 return bits_high[l1Id].bits[idx - bitct + L];
            });
            return skip_first_or_last_n_bits_and_count(v, bitId);
        }();

        auto right_l1 = (l1Id%2);
        auto right_l0 = (l0Id%2);
        auto symb0 = symb<<(bitct - L2);
        auto superblock = l0[l0Id/2][symb0-1];
        auto block      = l1[l1Id/2][symb0-1];
        auto r = [&]() {
            if (right_l0 && right_l1) {
                return superblock + block + count;
            } else if (right_l0) {
                return superblock + block - count;
            } else if (right_l1) {
                return superblock - block + count;
            } else {
                return superblock - block - count;
            }
        }();
        assert(r <= idx);
        return r;
    }

    auto prefix_rank_and_rank(uint64_t idx, uint64_t symb) const -> std::tuple<uint64_t, uint64_t> {
        return prefix_rank_and_rank_limit<bitct>(idx, symb);
    }

    template <size_t L2 = L>
    auto prefix_rank_and_rank_limit(uint64_t idx, uint64_t symb) const -> std::tuple<uint64_t, uint64_t> {
        assert(idx <= totalLength);
        assert(symb < size_t{1}<<L2);
        auto bitId = idx % (l1_bits_ct*2);
        auto l1Id = idx / l1_bits_ct;
        auto l0Id = idx / l0_bits_ct;
        assert(l1Id < bits_low.size());
        assert(l1Id/2 < l1.size());
        assert(l0Id/2 < l0.size());

        auto right_l1 = l1Id%2;
        auto right_l0 = l0Id%2;

        auto [superblock_pr, block_pr] = [&]() -> std::tuple<size_t, size_t> {
            if (symb == 0) return {0, 0};
            auto symb0 = symb<<(bitct - L2);
            return {l0[l0Id/2][symb0-1], l1[l1Id/2][symb0-1]};
        }();


        auto pr = [&]() -> size_t {
            if (symb == 0) return 0;
            if (symb == (size_t{1}<<L2)) return idx;


            auto count_pr = [&]() {
                auto v = mark_prefix_rank_cb<l1_bits_ct, L2>(symb, [&](size_t idx) -> auto const& {
                    idx += bitct - L2;
                    if (idx < bitct - L) return bits_low[l1Id].bits[idx];
                    else                 return bits_high[l1Id].bits[idx - bitct + L];
                });
                return skip_first_or_last_n_bits_and_count(v, bitId);
            }();
            if (right_l0 && right_l1) {
                return superblock_pr + block_pr + count_pr;
            } else if (right_l0) {
                return superblock_pr + block_pr - count_pr;
            } else if (right_l1) {
                return superblock_pr - block_pr + count_pr;
            } else {
                return superblock_pr - block_pr - count_pr;
            }
        }();

        auto count_r = [&]() {
            auto v = mark_prefix_rank_cb<l1_bits_ct, L2>(symb+1, [&](size_t idx) -> auto const& {
                idx += bitct - L2;
                if (idx < bitct - L) return bits_low[l1Id].bits[idx];
                else                 return bits_high[l1Id].bits[idx - bitct + L];
            });
            return skip_first_or_last_n_bits_and_count(v, bitId);
        }();


        auto symb1 = ((symb+1)<<(bitct-L2));
        auto superblock_r = l0[l0Id/2][symb1-1];
        auto block_r      = l1[l1Id/2][symb1-1];

        auto r = [&]() {
            if (right_l0 && right_l1) {
                return superblock_r + block_r + count_r;
            } else if (right_l0) {
                return superblock_r + block_r - count_r;
            } else if (right_l1) {
                return superblock_r - block_r + count_r;
            } else {
                return superblock_r - block_r - count_r;
            }
        }();

        assert(pr <= idx);
        assert(r <= idx);
        return {pr, r - pr};
    }

    auto all_ranks(uint64_t idx) const -> std::array<uint64_t, TSigma> {
        assert(idx <= totalLength);
        auto r = std::array<uint64_t, TSigma>{};
        for (size_t symb{0}; symb < TSigma; ++symb) {
            r[symb] = rank(idx, symb);
        }
        return r;
    }

    auto all_ranks_and_prefix_ranks(uint64_t idx) const -> std::tuple<std::array<uint64_t, TSigma>, std::array<uint64_t, TSigma>> {
        assert(idx <= totalLength);
        auto rs = all_ranks(idx);
        auto prs = std::array<uint64_t, TSigma>{};
        for (size_t i{1}; i < prs.size(); ++i) {
            prs[i] = prs[i-1] + rs[i-1];
        }
        return {rs, prs};
    }

    template <typename Archive>
    void serialize(this auto&& self, Archive& ar) {
        ar(self.l0, self.l1, self.bits_low, self.bits_high, self.totalLength);
    }
};

template <size_t Sigma> using PairedFlattenedBitvectors_b_64_4k   = PairedFlattenedBitvectors2L_b<Sigma, 2, 64, 4096>;
template <size_t Sigma> using PairedFlattenedBitvectors_b_128_4k  = PairedFlattenedBitvectors2L_b<Sigma, 2, 128, 4096>;
template <size_t Sigma> using PairedFlattenedBitvectors_b_256_4k  = PairedFlattenedBitvectors2L_b<Sigma, 2, 256, 4096>;
template <size_t Sigma> using PairedFlattenedBitvectors_b_512_4k  = PairedFlattenedBitvectors2L_b<Sigma, 2, 512, 4096>;
template <size_t Sigma> using PairedFlattenedBitvectors_b_1024_4k = PairedFlattenedBitvectors2L_b<Sigma, 2, 1024, 4096>;
template <size_t Sigma> using PairedFlattenedBitvectors_b_2048_4k = PairedFlattenedBitvectors2L_b<Sigma, 2, 2048, 4096>;

static_assert(checkStringKStep_c<PairedFlattenedBitvectors_b_64_4k>);
static_assert(checkStringKStep_c<PairedFlattenedBitvectors_b_128_4k>);
static_assert(checkStringKStep_c<PairedFlattenedBitvectors_b_256_4k>);
static_assert(checkStringKStep_c<PairedFlattenedBitvectors_b_512_4k>);
static_assert(checkStringKStep_c<PairedFlattenedBitvectors_b_1024_4k>);
static_assert(checkStringKStep_c<PairedFlattenedBitvectors_b_2048_4k>);

template <size_t Sigma> using PairedFlattenedBitvectors_b_64_64k   = PairedFlattenedBitvectors2L_b<Sigma, 2, 64, 65536>;
template <size_t Sigma> using PairedFlattenedBitvectors_b_128_64k  = PairedFlattenedBitvectors2L_b<Sigma, 2, 128, 65536>;
template <size_t Sigma> using PairedFlattenedBitvectors_b_256_64k  = PairedFlattenedBitvectors2L_b<Sigma, 2, 256, 65536>;
template <size_t Sigma> using PairedFlattenedBitvectors_b_512_64k  = PairedFlattenedBitvectors2L_b<Sigma, 2, 512, 65536>;
template <size_t Sigma> using PairedFlattenedBitvectors_b_1024_64k = PairedFlattenedBitvectors2L_b<Sigma, 2, 1024, 65536>;
template <size_t Sigma> using PairedFlattenedBitvectors_b_2048_64k = PairedFlattenedBitvectors2L_b<Sigma, 2, 2048, 65536>;
template <size_t Sigma> using PairedFlattenedBitvectors_b_4096_64k = PairedFlattenedBitvectors2L_b<Sigma, 2, 4096, 65536>;

static_assert(checkStringKStep_c<PairedFlattenedBitvectors_b_64_64k>);
static_assert(checkStringKStep_c<PairedFlattenedBitvectors_b_128_64k>);
static_assert(checkStringKStep_c<PairedFlattenedBitvectors_b_256_64k>);
static_assert(checkStringKStep_c<PairedFlattenedBitvectors_b_512_64k>);
static_assert(checkStringKStep_c<PairedFlattenedBitvectors_b_1024_64k>);
static_assert(checkStringKStep_c<PairedFlattenedBitvectors_b_2048_64k>);
static_assert(checkStringKStep_c<PairedFlattenedBitvectors_b_4096_64k>);

template <size_t Sigma> using PairedFlattenedBitvectors_b_64_64kUA   = PairedFlattenedBitvectors2L_b<Sigma, 2, 64, 65536, false>;
template <size_t Sigma> using PairedFlattenedBitvectors_b_128_64kUA  = PairedFlattenedBitvectors2L_b<Sigma, 2, 128, 65536, false>;
template <size_t Sigma> using PairedFlattenedBitvectors_b_256_64kUA  = PairedFlattenedBitvectors2L_b<Sigma, 2, 256, 65536, false>;
template <size_t Sigma> using PairedFlattenedBitvectors_b_512_64kUA  = PairedFlattenedBitvectors2L_b<Sigma, 2, 512, 65536, false>;
template <size_t Sigma> using PairedFlattenedBitvectors_b_1024_64kUA = PairedFlattenedBitvectors2L_b<Sigma, 2, 1024, 65536, false>;
template <size_t Sigma> using PairedFlattenedBitvectors_b_2048_64kUA = PairedFlattenedBitvectors2L_b<Sigma, 2, 2048, 65536, false>;
template <size_t Sigma> using PairedFlattenedBitvectors_b_4096_64kUA = PairedFlattenedBitvectors2L_b<Sigma, 2, 4096, 65536, false>;

static_assert(checkStringKStep_c<PairedFlattenedBitvectors_b_64_64kUA>);
static_assert(checkStringKStep_c<PairedFlattenedBitvectors_b_128_64kUA>);
static_assert(checkStringKStep_c<PairedFlattenedBitvectors_b_256_64kUA>);
static_assert(checkStringKStep_c<PairedFlattenedBitvectors_b_512_64kUA>);
static_assert(checkStringKStep_c<PairedFlattenedBitvectors_b_1024_64kUA>);
static_assert(checkStringKStep_c<PairedFlattenedBitvectors_b_2048_64kUA>);
static_assert(checkStringKStep_c<PairedFlattenedBitvectors_b_4096_64kUA>);

}
