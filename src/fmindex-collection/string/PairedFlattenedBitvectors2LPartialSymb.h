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

template <size_t TSigma, size_t K, size_t l1_bits_ct, size_t l0_bits_ct, bool Align=true>
struct PairedFlattenedBitvectors2LPartialSymb {
    static_assert(l1_bits_ct < l0_bits_ct, "first level must be smaller than second level");
    static_assert(l0_bits_ct-l1_bits_ct <= std::numeric_limits<uint16_t>::max(), "l0_bits_ct can only hold up to uint16_t bits");

    static constexpr size_t Sigma = TSigma;

    // number of full length bit vectors needed `2^bitct > TSigma`
    static constexpr auto bitct = std::bit_width(TSigma-1);
    // next full power of 2
    static constexpr auto bvct  = (1ull << bitct);

    static constexpr auto AlignVHigh = Align?alignAsValue(l1_bits_ct):alignof(std::array<std::bitset<l1_bits_ct>, K>);
    static constexpr auto AlignVLow  = Align?alignAsValue(l1_bits_ct):alignof(std::array<std::bitset<l1_bits_ct>, bitct<K?size_t{0}:(bitct-K)>);
    struct alignas(AlignVHigh) InBitsHigh {
        std::array<std::bitset<l1_bits_ct>, K> bits;

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
        std::array<std::bitset<l1_bits_ct>, bitct<K?size_t{0}:bitct-K> bits;

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

    PairedFlattenedBitvectors2LPartialSymb()
        : PairedFlattenedBitvectors2LPartialSymb{internal_tag{}, std::span<uint8_t const>{}}
    {}

    PairedFlattenedBitvectors2LPartialSymb(std::span<uint8_t const> _symbols)
        : PairedFlattenedBitvectors2LPartialSymb{internal_tag{}, _symbols}
    {}

    PairedFlattenedBitvectors2LPartialSymb(std::span<uint64_t const> _symbols)
        : PairedFlattenedBitvectors2LPartialSymb{internal_tag{}, _symbols}
    {}

    template <std::ranges::range range_t>
        requires std::convertible_to<std::ranges::range_value_t<range_t>, uint64_t>
    PairedFlattenedBitvectors2LPartialSymb(range_t&& _symbols)
        : PairedFlattenedBitvectors2LPartialSymb{internal_tag{}, _symbols}
    {}


private:
    struct internal_tag{};

    template <std::ranges::range range_t>
        requires std::convertible_to<std::ranges::range_value_t<range_t>, uint64_t>
    PairedFlattenedBitvectors2LPartialSymb(internal_tag, range_t&& _symbols) {

        if constexpr (requires() { _symbols.size(); }) {
            auto const _length = _symbols.size();
            bits_low.reserve(_length/l1_bits_ct + 2);
            bits_high.reserve(_length/l1_bits_ct + 2);
        }

        // fill all inbits
        for (auto c : _symbols) {
            auto bitId = totalLength % l1_bits_ct;
            for (size_t j{0}; j < bitct; ++j) {
                if (j < bitct - K) {
                    bits_low.back().bits[j][bitId] = (c>>j) & 1;
                } else {
                    bits_high.back().bits[j - bitct + K][bitId] = (c>>j) & 1;
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
                            if (idx < bitct - K) return b_low.bits[idx];
                            else                 return b_high.bits[idx-bitct+K];
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
                            if (idx < bitct - K) return b_low.bits[idx];
                            else                 return b_high.bits[idx-bitct+K];
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
        return symbol_limit<bitct>(idx);
    }

    template <size_t L = K>
    uint64_t symbol_limit(uint64_t idx) const {
        assert(idx < totalLength);
        auto bitId = idx % l1_bits_ct;
        auto l1Id  = idx / l1_bits_ct;
        assert(l1Id < bits_low.size());

        uint64_t symb{};
        for (uint64_t i{0}; i < L; ++i) {
            auto bit = [&](size_t i) {
                i += (bitct - L);
                if (i < bitct - K) return bits_low[l1Id].bits[i].test(bitId);
                else               return bits_high[l1Id].bits[i-bitct+K].test(bitId);
            }(i);
            symb = symb | (bit << i);
        }
        assert(symb < Sigma);
        return symb;
    }

    uint64_t rank(uint64_t idx, uint64_t symb) const {
        return rank_limit<bitct>(idx, symb);
    }

    template <size_t L = K>
    uint64_t rank_limit(uint64_t idx, uint64_t symb) const {
        assert(idx <= totalLength);
        assert(symb < (size_t{1}<<L));
        auto bitId = idx % (l1_bits_ct*2);
        auto l1Id = idx / l1_bits_ct;
        auto l0Id = idx / l0_bits_ct;
        assert(l1Id < bits_low.size());
        assert(l1Id/2 < l1.size());
        assert(l0Id/2 < l0.size());

        auto count = [&]() {
            auto v = mark_rank_cb<l1_bits_ct, L>(symb, [&](size_t idx) -> auto const& {
                idx += bitct - L;
                if (idx < bitct - K) return bits_low[l1Id].bits[idx];
                else                 return bits_high[l1Id].bits[idx-bitct+K];
            });
            return skip_first_or_last_n_bits_and_count(v, bitId);
        }();

        auto symb1 = ((symb+1)<<(bitct-L))-1;
        auto right_l1 = (l1Id%2);
        auto right_l0 = (l0Id%2);
        auto [superblock, block] = [&]() -> std::tuple<size_t, size_t> {
            if (symb == 0) return {l0[l0Id/2][symb1], l1[l1Id/2][symb1]};
            auto symb0 = (symb<<(bitct-L))-1;
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

    template <size_t L = K>
    uint64_t prefix_rank_limit(uint64_t idx, uint64_t symb) const {
        if (symb == 0) return 0;
        if (symb == (size_t{1}<<L)) return idx;
        assert(idx <= totalLength);
        assert(symb <= (size_t{1}<<L));
        auto bitId = idx % (l1_bits_ct*2);
        auto l1Id = idx / l1_bits_ct;
        auto l0Id = idx / l0_bits_ct;
        assert(l1Id < bits_low.size());
        assert(l1Id/2 < l1.size());
        assert(l0Id/2 < l0.size());

        auto count = [&]() {
            auto v = mark_prefix_rank_cb<l1_bits_ct, L>(symb, [&](size_t idx) -> auto const& {
                idx += bitct - L;
                if (idx < bitct - K) return bits_low[l1Id].bits[idx];
                else                 return bits_high[l1Id].bits[idx - bitct + K];
            });
            return skip_first_or_last_n_bits_and_count(v, bitId);
        }();

        auto right_l1 = (l1Id%2);
        auto right_l0 = (l0Id%2);
        auto symb0 = symb<<(bitct - L);
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

    template <size_t L = K>
    auto prefix_rank_and_rank_limit(uint64_t idx, uint64_t symb) const -> std::tuple<uint64_t, uint64_t> {
        assert(idx <= totalLength);
        assert(symb < size_t{1}<<L);
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
            auto symb0 = symb<<(bitct - L);
            return {l0[l0Id/2][symb0-1], l1[l1Id/2][symb0-1]};
        }();


        auto pr = [&]() -> size_t {
            if (symb == 0) return 0;
            if (symb == (size_t{1}<<L)) return idx;


            auto count_pr = [&]() {
                auto v = mark_prefix_rank_cb<l1_bits_ct, L>(symb, [&](size_t idx) -> auto const& {
                    idx += bitct - L;
                    if (idx < bitct - K) return bits_low[l1Id].bits[idx];
                    else                 return bits_high[l1Id].bits[idx - bitct + K];
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
            auto v = mark_prefix_rank_cb<l1_bits_ct, L>(symb+1, [&](size_t idx) -> auto const& {
                idx += bitct - L;
                if (idx < bitct - K) return bits_low[l1Id].bits[idx];
                else                 return bits_high[l1Id].bits[idx - bitct + K];
            });
            return skip_first_or_last_n_bits_and_count(v, bitId);
        }();


        auto symb1 = ((symb+1)<<(bitct-L));
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

    /* This is a special variant of the `all_ranks_and_prefix_ranks` function.
     *
     * Differences:
     *  - it looks up the values for two positions at the same time.
     *  - it uses a callback function
     *  - it does not report (prefix) ranks for symbols that have 0-values
     */
    template <size_t L>
    void all_ranks_dual_limit(size_t idx1, size_t idx2, auto const& cb) const {
        assert(idx1 <= totalLength);
        assert(idx2 <= totalLength);

        // takes an index and splits it into the indices required for partial-paired-blocks
        auto split = [&](size_t idx) -> std::tuple<size_t, size_t, size_t, int64_t, int64_t> {
            auto bitId = idx % (l1_bits_ct*2);
            auto l1Id = idx / l1_bits_ct;
            auto l0Id = idx / l0_bits_ct;
            assert(l1Id < bits_low.size());
            assert(l1Id/2 < l1.size());
            assert(l0Id/2 < l0.size());
            int64_t right_l1 = l1Id%2;
            int64_t right_l0 = l0Id%2;
            assert(bitId <= l1_bits_ct*2);
            return {bitId, l1Id, l0Id, right_l1, right_l0};
        };

        // compute indices for idx1 (_lb) and idx2 (_rb)
        auto [bitId_lb, l1Id_lb, l0Id_lb, right_l1_lb, right_l0_lb] = split(idx1);
        auto [bitId_rb, l1Id_rb, l0Id_rb, right_l1_rb, right_l0_rb] = split(idx2);

        // compute blocks
        auto const& l0_lb = l0[l0Id_lb/2];
        auto const& l0_rb = l0[l0Id_rb/2];
        auto const& l1_lb = l1[l1Id_lb/2];
        auto const& l1_rb = l1[l1Id_rb/2];

        auto bits_lb = [&](size_t idx) {
            if (idx < bitct - K) return bits_low[l1Id_lb].bits[idx];
            else                 return bits_high[l1Id_lb].bits[idx - bitct + K];
        };
        auto bits_rb = [&](size_t idx) {
            if (idx < bitct - K) return bits_low[l1Id_rb].bits[idx];
            else                 return bits_high[l1Id_rb].bits[idx - bitct + K];
        };
//        auto const& bits_lb = bits[l1Id_lb].bits;
//        auto const& bits_rb = bits[l1Id_rb].bits;


        // computes prefix rank for symbol (symb)
        auto count_pr = [&](auto const& b1, auto const& b2, size_t symb) -> std::tuple<size_t, size_t> {
            auto count_lb = skip_first_or_last_n_bits_and_count(b1, bitId_lb);
            auto count_rb = skip_first_or_last_n_bits_and_count(b2, bitId_rb);

            auto pr1 = [&]() {
                auto const& superblock = (symb>0)?l0_lb[symb-1]:size_t{0};
                auto const& block      = (symb>0)?l1_lb[symb-1]:size_t{0};
                auto const& count      = count_lb;
                auto const& right_l0   = right_l0_lb;
                auto const& right_l1   = right_l1_lb;

                if (right_l0 && right_l1) return superblock + block + count;
                else if (right_l0)        return superblock + block - count;
                else if (right_l1)        return superblock - block + count;
                else                      return superblock - block - count;
            }();

            auto pr2 = [&]() {
                auto const& superblock = (symb>0)?l0_rb[symb-1]:size_t{0};
                auto const& block      = (symb>0)?l1_rb[symb-1]:size_t{0};
                auto const& count      = count_rb;
                auto const& right_l0   = right_l0_rb;
                auto const& right_l1   = right_l1_rb;

                if (right_l0 && right_l1) return superblock + block + count;
                else if (right_l0)        return superblock + block - count;
                else if (right_l1)        return superblock - block + count;
                else                      return superblock - block - count;
            }();

            assert(pr1 <= idx1);
            assert(pr2 <= idx2);

            return {pr1, pr2};
        };

        using word_of_bits = std::bitset<l1_bits_ct>;
        auto This = this; //!TODO !WORKAROUND otherwise crashes gcc
        (void)This;

        /** Heart of the algorithm
         *
         * A recursive algorithm checks if descending to left and/or right is required.
         * Aborts if any of the ranges have no characters. This accelerates computation on large alphabets.
         */
        auto rec = [&](this auto&& self, word_of_bits const& l_b1, word_of_bits const& r_b1, word_of_bits const& l_b2, word_of_bits const& r_b2, int level, size_t pr1_S, size_t pr2_S, size_t pr1_E, size_t pr2_E, size_t symb) -> void {
            auto lsymb = symb | (1 << level);
            auto b1 = l_b1 | (r_b1 & ~bits_lb(level));
            auto b2 = l_b2 | (r_b2 & ~bits_rb(level));

            auto [pr1, pr2] = count_pr(b1, b2, lsymb);

            assert(This->prefix_rank(idx1, lsymb) == pr1);
            assert(This->prefix_rank(idx2, lsymb) == pr2);

            auto rank_l_lb = pr1-pr1_S;
            auto rank_l_rb = pr2-pr2_S;

            auto rank_r_lb = pr1_E-pr1;
            auto rank_r_rb = pr2_E-pr2;
            if (level == bitct-L) {
                assert(This->prefix_rank(idx1, symb) == pr1_S);
                assert(This->prefix_rank(idx2, symb) == pr2_S);
                cb( symb >> (bitct-L), rank_l_lb, rank_l_rb, pr1_S, pr2_S);
                cb(lsymb >> (bitct-L), rank_r_lb, rank_r_rb, pr1, pr2);
            } else {
                auto countLeft  = rank_l_rb - rank_l_lb;
                auto countRight = rank_r_rb - rank_r_lb;

                // recursion left
                if (countLeft > 0) {
                    self(l_b1, b1, l_b2, b2, level-1, pr1_S, pr2_S, pr1, pr2, symb);
                }
                // recursion right
                if (countRight > 0) {
                    self(b1, r_b1, b2, r_b2, level-1, pr1, pr2, pr1_E, pr2_E, lsymb);
                }
            }
        };
        rec (
            mask_positive_or_negative<l1_bits_ct>[1], mask_positive_or_negative<l1_bits_ct>[0],
            mask_positive_or_negative<l1_bits_ct>[1], mask_positive_or_negative<l1_bits_ct>[0],
            /*.level=*/      bitct-1,
            0, 0, idx1, idx2,
            /*.symb=*/       0
        );
    }

    template <typename Archive>
    void serialize(this auto&& self, Archive& ar) {
        ar(self.l0, self.l1, self.bits_low, self.bits_high, self.totalLength);
    }
};

template <size_t Sigma> using PairedFlattenedBitvectorsPartialSymb_64_4k   = PairedFlattenedBitvectors2LPartialSymb<Sigma, 2, 64, 4096>;
template <size_t Sigma> using PairedFlattenedBitvectorsPartialSymb_128_4k  = PairedFlattenedBitvectors2LPartialSymb<Sigma, 2, 128, 4096>;
template <size_t Sigma> using PairedFlattenedBitvectorsPartialSymb_256_4k  = PairedFlattenedBitvectors2LPartialSymb<Sigma, 2, 256, 4096>;
template <size_t Sigma> using PairedFlattenedBitvectorsPartialSymb_512_4k  = PairedFlattenedBitvectors2LPartialSymb<Sigma, 2, 512, 4096>;
template <size_t Sigma> using PairedFlattenedBitvectorsPartialSymb_1024_4k = PairedFlattenedBitvectors2LPartialSymb<Sigma, 2, 1024, 4096>;
template <size_t Sigma> using PairedFlattenedBitvectorsPartialSymb_2048_4k = PairedFlattenedBitvectors2LPartialSymb<Sigma, 2, 2048, 4096>;

static_assert(checkStringKStep_c<PairedFlattenedBitvectorsPartialSymb_64_4k>);
static_assert(checkStringKStep_c<PairedFlattenedBitvectorsPartialSymb_128_4k>);
static_assert(checkStringKStep_c<PairedFlattenedBitvectorsPartialSymb_256_4k>);
static_assert(checkStringKStep_c<PairedFlattenedBitvectorsPartialSymb_512_4k>);
static_assert(checkStringKStep_c<PairedFlattenedBitvectorsPartialSymb_1024_4k>);
static_assert(checkStringKStep_c<PairedFlattenedBitvectorsPartialSymb_2048_4k>);

template <size_t Sigma> using PairedFlattenedBitvectorsPartialSymb_64_64k   = PairedFlattenedBitvectors2LPartialSymb<Sigma, 2, 64, 65536>;
template <size_t Sigma> using PairedFlattenedBitvectorsPartialSymb_128_64k  = PairedFlattenedBitvectors2LPartialSymb<Sigma, 2, 128, 65536>;
template <size_t Sigma> using PairedFlattenedBitvectorsPartialSymb_256_64k  = PairedFlattenedBitvectors2LPartialSymb<Sigma, 2, 256, 65536>;
template <size_t Sigma> using PairedFlattenedBitvectorsPartialSymb_512_64k  = PairedFlattenedBitvectors2LPartialSymb<Sigma, 2, 512, 65536>;
template <size_t Sigma> using PairedFlattenedBitvectorsPartialSymb_1024_64k = PairedFlattenedBitvectors2LPartialSymb<Sigma, 2, 1024, 65536>;
template <size_t Sigma> using PairedFlattenedBitvectorsPartialSymb_2048_64k = PairedFlattenedBitvectors2LPartialSymb<Sigma, 2, 2048, 65536>;
template <size_t Sigma> using PairedFlattenedBitvectorsPartialSymb_4096_64k = PairedFlattenedBitvectors2LPartialSymb<Sigma, 2, 4096, 65536>;

static_assert(checkStringKStep_c<PairedFlattenedBitvectorsPartialSymb_64_64k>);
static_assert(checkStringKStep_c<PairedFlattenedBitvectorsPartialSymb_128_64k>);
static_assert(checkStringKStep_c<PairedFlattenedBitvectorsPartialSymb_256_64k>);
static_assert(checkStringKStep_c<PairedFlattenedBitvectorsPartialSymb_512_64k>);
static_assert(checkStringKStep_c<PairedFlattenedBitvectorsPartialSymb_1024_64k>);
static_assert(checkStringKStep_c<PairedFlattenedBitvectorsPartialSymb_2048_64k>);
static_assert(checkStringKStep_c<PairedFlattenedBitvectorsPartialSymb_4096_64k>);

template <size_t Sigma> using PairedFlattenedBitvectorsPartialSymb_64_64kUA   = PairedFlattenedBitvectors2LPartialSymb<Sigma, 2, 64, 65536, false>;
template <size_t Sigma> using PairedFlattenedBitvectorsPartialSymb_128_64kUA  = PairedFlattenedBitvectors2LPartialSymb<Sigma, 2, 128, 65536, false>;
template <size_t Sigma> using PairedFlattenedBitvectorsPartialSymb_256_64kUA  = PairedFlattenedBitvectors2LPartialSymb<Sigma, 2, 256, 65536, false>;
template <size_t Sigma> using PairedFlattenedBitvectorsPartialSymb_512_64kUA  = PairedFlattenedBitvectors2LPartialSymb<Sigma, 2, 512, 65536, false>;
template <size_t Sigma> using PairedFlattenedBitvectorsPartialSymb_1024_64kUA = PairedFlattenedBitvectors2LPartialSymb<Sigma, 2, 1024, 65536, false>;
template <size_t Sigma> using PairedFlattenedBitvectorsPartialSymb_2048_64kUA = PairedFlattenedBitvectors2LPartialSymb<Sigma, 2, 2048, 65536, false>;
template <size_t Sigma> using PairedFlattenedBitvectorsPartialSymb_4096_64kUA = PairedFlattenedBitvectors2LPartialSymb<Sigma, 2, 4096, 65536, false>;

static_assert(checkStringKStep_c<PairedFlattenedBitvectorsPartialSymb_64_64kUA>);
static_assert(checkStringKStep_c<PairedFlattenedBitvectorsPartialSymb_128_64kUA>);
static_assert(checkStringKStep_c<PairedFlattenedBitvectorsPartialSymb_256_64kUA>);
static_assert(checkStringKStep_c<PairedFlattenedBitvectorsPartialSymb_512_64kUA>);
static_assert(checkStringKStep_c<PairedFlattenedBitvectorsPartialSymb_1024_64kUA>);
static_assert(checkStringKStep_c<PairedFlattenedBitvectorsPartialSymb_2048_64kUA>);
static_assert(checkStringKStep_c<PairedFlattenedBitvectorsPartialSymb_4096_64kUA>);

}
