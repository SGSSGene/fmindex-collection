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


template <size_t TSigma, size_t l1_bits_ct, size_t l0_bits_ct, bool Align=true>
struct PairedFlattenedBitvectors2L {
    static_assert(l1_bits_ct < l0_bits_ct, "first level must be smaller than second level");
    static_assert(l0_bits_ct-l1_bits_ct <= std::numeric_limits<uint16_t>::max(), "l0_bits_ct can only hold up to uint16_t bits");

    static constexpr size_t Sigma = TSigma;

    // number of full length bit vectors needed `2^bitct > TSigma`
    static constexpr auto bitct = std::bit_width(TSigma-1);
    // next full power of 2
    static constexpr auto bvct  = (1ull << bitct);

    static constexpr auto AlignV = Align?alignAsValue(l1_bits_ct):alignof(std::array<std::bitset<l1_bits_ct>, bitct>);
    struct alignas(AlignV) InBits {
        std::array<std::bitset<l1_bits_ct>, bitct> bits;

        uint64_t symbol(uint64_t idx) const {
            assert(idx < l1_bits_ct);
            uint64_t symb{};
            for (uint64_t i{bitct}; i > 0; --i) {
                auto b = bits[i-1].test(idx);
                symb = (symb<<1) | b;
            }
            assert(symb < Sigma);
            return symb;
        }

        uint64_t rank(uint64_t idx, uint64_t symb) const {
            assert(idx <= l1_bits_ct*2);
            assert(symb < Sigma);
            auto v = neprv8_detail::rank(std::span{bits}, symb);
            return skip_first_or_last_n_bits_and_count(v, idx);
        }

        template <size_t L>
        uint64_t rank_limit(uint64_t idx, uint64_t symb) const {
            assert(idx <= l1_bits_ct*2);
            assert(symb < (size_t{1}<<L));
            auto bits_span = std::span{bits}.template last<L>();
            auto v = neprv8_detail::rank(bits_span, symb);
            return skip_first_or_last_n_bits_and_count(v, idx);
        }

        uint64_t prefix_rank(uint64_t idx, uint64_t symb) const {
            assert(idx <= l1_bits_ct*2);
            assert(symb <= Sigma);
            auto v = neprv8_detail::prefix_rank(std::span{bits}, symb);
            return skip_first_or_last_n_bits_and_count(v, idx);
        }

        template <size_t L>
        uint64_t prefix_rank_limit(uint64_t idx, uint64_t symb) const {
            assert(idx <= l1_bits_ct*2);
            assert(symb < (size_t{1}<<L));
            auto bits_span = std::span{bits}.template last<L>();
            auto v = neprv8_detail::prefix_rank(bits_span, symb);
            return skip_first_or_last_n_bits_and_count(v, idx);
        }

        auto all_ranks(uint64_t idx) const -> std::array<uint64_t, TSigma> {
            assert(idx <= l1_bits_ct*2);

            auto vs = neprv8_detail::rank_all<(1ull<<bitct)>(std::span{bits});
            auto v = std::array<uint64_t, TSigma>{};
            static_assert(v.size() <= vs.size());
            for (size_t i{0}; i < v.size(); ++i) {
                auto count = skip_first_or_last_n_bits_and_count(vs[i], idx);
                v[i] = count;
            }
            return v;
        }
        void setSymbol(size_t i, uint64_t symb) {
            assert(i < l1_bits_ct);
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

    mmser::vector<InBits> bits{{}};
    mmser::vector<BlockL1> l1{{}};
    mmser::vector<BlockL0> l0{{}};

/*    using BlockL1Small = std::array<uint16_t, 4>;
    using BlockL0Small = std::array<uint64_t, 4>;
    mmser::vector<BlockL1Small> l1_small{{}};
    mmser::vector<BlockL0Small> l0_small{{}};*/


    size_t totalLength{};

    PairedFlattenedBitvectors2L()
        : PairedFlattenedBitvectors2L{internal_tag{}, std::span<uint8_t const>{}}
    {}

    PairedFlattenedBitvectors2L(std::span<uint8_t const> _symbols)
        : PairedFlattenedBitvectors2L{internal_tag{}, _symbols}
    {}

    PairedFlattenedBitvectors2L(std::span<uint64_t const> _symbols)
        : PairedFlattenedBitvectors2L{internal_tag{}, _symbols}
    {}

    template <std::ranges::range range_t>
        requires std::convertible_to<std::ranges::range_value_t<range_t>, uint64_t>
    PairedFlattenedBitvectors2L(range_t&& _symbols)
        : PairedFlattenedBitvectors2L{internal_tag{}, _symbols}
    {}


private:
    struct internal_tag{};

    template <std::ranges::range range_t>
        requires std::convertible_to<std::ranges::range_value_t<range_t>, uint64_t>
    PairedFlattenedBitvectors2L(internal_tag, range_t&& _symbols) {

        if constexpr (requires() { _symbols.size(); }) {
            auto const _length = _symbols.size();
            bits.reserve(_length/l1_bits_ct + 2);
        }

        // fill all inbits
        for (auto c : _symbols) {
            auto bitId = totalLength % l1_bits_ct;
            bits.back().setSymbol(bitId, c);

            totalLength += 1;
            if (totalLength % l1_bits_ct == 0) { // next bit will require a new in-bits block
                bits.emplace_back();
            }
        }

        // fill l0/l1 structure
        {
            size_t l0BlockCt = (totalLength / (l0_bits_ct*2)) + 1;
            size_t l1BlockCt = l0BlockCt * (l0_bits_ct / l1_bits_ct);
            size_t inbitsCt  = l0BlockCt * ((l0_bits_ct*2) / l1_bits_ct);

            l0.resize(l0BlockCt);
            l1.resize(l1BlockCt);
            bits.resize(inbitsCt);

            constexpr size_t b1 = l1_bits_ct;
            constexpr size_t b0 = l0_bits_ct;

            constexpr size_t l1_block_ct = b0 / b1;

            BlockL0 l0_acc{};
            // walk through all superblocks
            for (size_t l0I{0}; l0I < l0BlockCt; ++l0I) {
                // walk left to right and set l1 values (as if they are the begining of a superblock)

                // left part
                BlockL0 acc{};
                for (size_t i{0}; i < l1_block_ct; ++i) {
                    auto idx = l0I*l1_block_ct*2 + i;
                    auto& b = bits[idx];
                    auto counts = b.all_ranks(0);
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
                    auto& b = bits[idx];
                    auto counts = b.all_ranks(0);

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

/*            if constexpr (Sigma == 16) {
                for (auto& v : l0.owningBuffer) {
                    BlockL0Small a;
                    a[0] = v[0*4+0];
                    a[1] = v[1*4+0];
                    a[2] = v[2*4+0];
                    a[3] = v[3*4+0];
                    l0_small.push_back(a);
                }
                for (auto& v : l1.owningBuffer) {
                    BlockL1Small a;
                    a[0] = v[0*4+0];
                    a[1] = v[1*4+0];
                    a[2] = v[2*4+0];
                    a[3] = v[3*4+0];
                    l1_small.push_back(a);
                }
            }*/



/*            if constexpr (Sigma == 16) {
                for (auto& v : l0.owningBuffer) {
                    BlockL0 arr;
                    arr[ 0] = v[0*4+0];
                    arr[ 1] = v[1*4+0];
                    arr[ 2] = v[2*4+0];
                    arr[ 3] = v[3*4+0];
                    arr[ 4] = v[0*4+1];
                    arr[ 5] = v[1*4+1];
                    arr[ 6] = v[2*4+1];
                    arr[ 7] = v[3*4+1];
                    arr[ 8] = v[0*4+2];
                    arr[ 9] = v[1*4+2];
                    arr[10] = v[2*4+2];
                    arr[11] = v[3*4+2];
                    arr[12] = v[0*4+3];
                    arr[13] = v[1*4+3];
                    arr[14] = v[2*4+3];
                    arr[15] = v[3*4+3];
                    v = arr;
                }
                for (auto& v : l1.owningBuffer) {
                    BlockL1 arr;
                    arr[ 0] = v[0*4+0];
                    arr[ 1] = v[1*4+0];
                    arr[ 2] = v[2*4+0];
                    arr[ 3] = v[3*4+0];
                    arr[ 4] = v[0*4+1];
                    arr[ 5] = v[1*4+1];
                    arr[ 6] = v[2*4+1];
                    arr[ 7] = v[3*4+1];
                    arr[ 8] = v[0*4+2];
                    arr[ 9] = v[1*4+2];
                    arr[10] = v[2*4+2];
                    arr[11] = v[3*4+2];
                    arr[12] = v[0*4+3];
                    arr[13] = v[1*4+3];
                    arr[14] = v[2*4+3];
                    arr[15] = v[3*4+3];
                    v = arr;
                }
            }*/
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
        assert(l1Id < bits.size());
        auto symb = bits[l1Id].symbol(bitId);
        assert(symb < Sigma);
        return symb;
    }

    template <size_t L>
    uint64_t symbol_limit(uint64_t idx) const {
        return (symbol(idx) >> (bitct-L)); //!TODO
    }

/*    template <size_t L>
    uint64_t symbol_tail(uint64_t idx) const {
        return symbol(idx) & ((size_t{1} << (L)) - 1);
    }*/


    uint64_t rank(uint64_t idx, uint64_t symb) const {
        assert(idx <= totalLength);
        assert(symb < Sigma);
        auto bitId = idx % (l1_bits_ct*2);
        auto l1Id = idx / l1_bits_ct;
        auto l0Id = idx / l0_bits_ct;
        assert(l1Id < bits.size());
        assert(l1Id/2 < l1.size());
        assert(l0Id/2 < l0.size());

        auto count      = bits[l1Id].rank(bitId, symb);
        auto right_l1   = (l1Id%2);
        auto right_l0   = (l0Id%2);
        auto [superblock, block] = [&]() -> std::tuple<size_t, size_t> {
            if (symb == 0) return {l0[l0Id/2][0], l1[l1Id/2][0]};
            return {
                (l0[l0Id/2][symb] - l0[l0Id/2][symb-1]),
                (l1[l1Id/2][symb] - l1[l1Id/2][symb-1])
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

    template <size_t L>
    uint64_t rank_limit(uint64_t idx, uint64_t symb) const {
        assert(idx <= totalLength);
        assert(symb < (size_t{1}<<L));
        auto bitId = idx % (l1_bits_ct*2);
        auto l1Id = idx / l1_bits_ct;
        auto l0Id = idx / l0_bits_ct;
        assert(l1Id < bits.size());
        assert(l1Id/2 < l1.size());
        assert(l0Id/2 < l0.size());

        auto count = bits[l1Id].template rank_limit<L>(bitId, symb);
/*        auto count1 = bits[l1Id].template prefix_rank_limit<L>(bitId, symb);
        auto count2 = bits[l1Id].template prefix_rank_limit<L>(bitId, symb+1);
        auto count = count2 - count1;*/


        if constexpr (bitct == 4 && L == 2 && false) {
/*             auto symb0 = symb;
             auto symb1 = symb+1;
             (void)symb1;

            auto right_l1 = (l1Id%2);
            auto right_l0 = (l0Id%2);
            auto superblock = (l0_small[l0Id/2][symb1] - l0_small[l0Id/2][symb0]);
            auto block      = (l1_small[l1Id/2][symb1] - l1_small[l1Id/2][symb0]);
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
            return r;*/
        } else {
//             auto symb0 = symb;
//             auto symb1 = symb+1;
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
    }


    uint64_t prefix_rank(uint64_t idx, uint64_t symb) const {
        if (symb == 0) return 0;
        if (symb == Sigma) return idx;
        assert(idx <= totalLength);
        assert(symb <= Sigma);
        auto bitId = idx % (l1_bits_ct*2);
        auto l1Id = idx / l1_bits_ct;
        auto l0Id = idx / l0_bits_ct;
        assert(l1Id < bits.size());
        assert(l1Id/2 < l1.size());
        assert(l0Id/2 < l0.size());

        auto count = bits[l1Id].prefix_rank(bitId, symb);

        auto right_l1 = (l1Id%2);
        auto right_l0 = (l0Id%2);
        auto superblock = l0[l0Id/2][symb-1];
        auto block      = l1[l1Id/2][symb-1];
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

    template <size_t L>
    uint64_t prefix_rank_limit(uint64_t idx, uint64_t symb) const {
        if (symb == 0) return 0;
        if (symb == Sigma) return idx;
        assert(idx <= totalLength);
        assert(symb <= (size_t{1}<<L));
        auto bitId = idx % (l1_bits_ct*2);
        auto l1Id = idx / l1_bits_ct;
        auto l0Id = idx / l0_bits_ct;
        assert(l1Id < bits.size());
        assert(l1Id/2 < l1.size());
        assert(l0Id/2 < l0.size());

        auto count = bits[l1Id].template prefix_rank_limit<L>(bitId, symb);

/*        if constexpr (bitct == 4 && L == 2 && true) {
            auto right_l1 = (l1Id%2);
            auto right_l0 = (l0Id%2);
            auto superblock = l0_small[l0Id/2][symb-1];
            auto block      = l1_small[l1Id/2][symb-1];
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

        } else {*/
            auto symb0 = symb<<(bitct - L);

            auto right_l1 = (l1Id%2);
            auto right_l0 = (l0Id%2);
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
//        }
    }

    auto prefix_rank_and_rank(uint64_t idx, uint64_t symb) const -> std::tuple<size_t, size_t> {
/*        auto pr0 = prefix_rank(idx, symb);
        auto pr1 = prefix_rank(idx, symb+1);
        return {pr0, pr1 - pr0};*/

//        return {prefix_rank(idx, symb), rank(idx, symb)};
        assert(idx <= totalLength);
        assert(symb < Sigma);
        auto bitId = idx % (l1_bits_ct*2);
        auto l1Id = idx / l1_bits_ct;
        auto l0Id = idx / l0_bits_ct;
        assert(l1Id < bits.size());
        assert(l1Id/2 < l1.size());
        assert(l0Id/2 < l0.size());

        __builtin_prefetch(reinterpret_cast<void const*>(&l0[l0Id/2][symb]), 0, 0);
        __builtin_prefetch(reinterpret_cast<void const*>(&l1[l1Id/2][symb]), 0, 0);
        __builtin_prefetch(reinterpret_cast<void const*>(&bits[l1Id]), 0, 0);

/*        for (auto& b : bits[l1Id].bits) {
            __builtin_prefetch(reinterpret_cast<void const*>(&b), 0, 0);
        }*/

        auto pr0 = prefix_rank(idx, symb);
        auto pr1 = prefix_rank(idx, symb+1);
        return {pr0, pr1 - pr0};

        auto right_l1 = l1Id%2;
        auto right_l0 = l0Id%2;


        auto [superblock_pr, block_pr] = [&]() -> std::tuple<size_t, size_t> {
            if (symb == 0) return {0, 0};
            auto symb0 = symb - 1;
            return {l0[l0Id/2][symb0], l1[l1Id/2][symb0]};
        }();


        auto pr = [&]() -> size_t {
            if (symb == 0) return 0;

            auto count_pr = bits[l1Id].prefix_rank(bitId, symb);
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

//        auto count_r = bits[l1Id].rank(bitId, symb);
        auto count_r = bits[l1Id].prefix_rank(bitId, symb+1);

        auto symb1 = symb;
        auto superblock_r = l0[l0Id/2][symb1];
        auto block_r      = l1[l1Id/2][symb1];

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

    template <size_t L>
    auto prefix_rank_and_rank_limit(uint64_t idx, uint64_t symb) const -> std::tuple<size_t, size_t> {
        {
            auto l1Id = idx / l1_bits_ct;
            auto l0Id = idx / l0_bits_ct;
            __builtin_prefetch(reinterpret_cast<void const*>(&l0[l0Id/2][symb<<(bitct-L)]), 0, 0);
            __builtin_prefetch(reinterpret_cast<void const*>(&l1[l1Id/2][symb<<(bitct-L)]), 0, 0);
            __builtin_prefetch(reinterpret_cast<void const*>(&bits[l1Id].bits.back()), 0, 0);
/*            for (size_t i{L}; i < bitct; ++i) {
                auto& b = bits[l1Id].bits[i];
                __builtin_prefetch(reinterpret_cast<void const*>(&b), 0, 0);
            }*/
        }

        auto pr0 = prefix_rank_limit<L>(idx, symb);
        auto pr1 = prefix_rank_limit<L>(idx, symb+1);
        return {pr0, pr1 - pr0};

        return {prefix_rank_limit<L>(idx, symb), rank_limit<L>(idx, symb)};
        assert(idx <= totalLength);
        assert(symb < (size_t{1}<<L));
        auto bitId = idx % (l1_bits_ct*2);
        auto l1Id = idx / l1_bits_ct;
        auto l0Id = idx / l0_bits_ct;
        assert(l1Id < bits.size());
        assert(l1Id/2 < l1.size());
        assert(l0Id/2 < l0.size());

/*        auto count_r  = bits[l1Id].template rank_limit<L>(bitId, symb);
        auto count_pr = bits[l1Id].template prefix_rank_limit<L>(bitId, symb);*/

/*        auto count_pr = bits[l1Id].prefix_rank(bitId, symb);
        auto count_r  = bits[l1Id].prefix_rank(bitId, symb+(size_t{1ull}<<L)) - count_pr;*/

        auto count_r  = bits[l1Id].template rank_limit<L>(bitId, symb);
        auto count_pr = size_t{0};
//        auto const step = (size_t{1}<<L);
        for (size_t i{0}; i < symb; i += 1) {
            count_pr += bits[l1Id].template rank_limit<L>(bitId, symb);
        }


        auto right_l1 = l1Id%2;
        auto right_l0 = l0Id%2;

        auto symb1 = ((symb+1)<<(bitct-L))-1;
        auto [superblock_pr, block_pr] = [&]() -> std::tuple<size_t, size_t> {
            if (symb == 0) return {0, 0};
            auto symb0 = (symb<<(bitct-L))-1;
            return {l0[l0Id/2][symb0], l1[l1Id/2][symb0]};
        }();
        auto [superblock_r, block_r] = [&]() -> std::tuple<size_t, size_t> {
            return {l0[l0Id/2][symb1] - superblock_pr, l1[l1Id/2][symb1] - block_pr};
        }();

        auto [pr, r] = [&]() -> std::tuple<size_t, size_t> {
            if (right_l0 && right_l1) {
                return {
                    superblock_pr + block_pr + count_pr,
                    superblock_r  + block_r  + count_r
                };
            } else if (right_l0) {
                return {
                    superblock_pr + block_pr - count_pr,
                    superblock_r  + block_r  - count_r
                };
            } else if (right_l1) {
                return {
                    superblock_pr - block_pr + count_pr,
                    superblock_r  - block_r  + count_r
                };
            } else {
                return {
                    superblock_pr - block_pr - count_pr,
                    superblock_r  - block_r  - count_r
                };
            }
        }();

        assert(pr <= idx);
        assert(r <= idx);
        return {pr, r};
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
    void all_ranks_dual(size_t idx1, size_t idx2, auto const& cb) const {
        assert(idx1 <= totalLength);
        assert(idx2 <= totalLength);

        /*
         * The below is a "optimized version" of the following easy code
            size_t prevPrefixRank1{};
            size_t prevPrefixRank2{};
            for (size_t symb{0}; symb < TSigma; ++symb) {
                auto pr1 = prefix_rank(idx1, symb+1);
                auto pr2 = prefix_rank(idx2, symb+1);
                cb(symb, pr1-prevPrefixRank1, pr2-prevPrefixRank2, prevPrefixRank1, prevPrefixRank2);
                prevPrefixRank1 = pr1;
                prevPrefixRank2 = pr2;
            }
        */

        // takes an index and splits it into the indices required for partial-paired-blocks
        auto split = [&](size_t idx) -> std::tuple<size_t, size_t, size_t, size_t, size_t> {
            auto bitId = idx % (l1_bits_ct*2);
            auto l1Id = idx / l1_bits_ct;
            auto l0Id = idx / l0_bits_ct;
            assert(l1Id < bits.size());
            assert(l1Id/2 < l1.size());
            assert(l0Id/2 < l0.size());
            auto right_l1 = l1Id%2;
            auto right_l0 = l0Id%2;
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
        auto const& bits_lb = bits[l1Id_lb].bits;
        auto const& bits_rb = bits[l1Id_rb].bits;

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



        /*

        0      ~b3 & ~b2 & ~b1 | ~b0
        1      ~b3 & ~b2 & ~b1
        2      ~b3 & ~b2 & (~b1 | ~b0)
        3
        4
        5
        6
        7


        0
        1         ~b1 & ~b0
        2    ~b1
        3         ~b1 | ~b0



  000   0
  001   1                      ~b2 & ~b1 & ~b0
  010   2        ~b2 & ~b1
  011   3                      (~b2 & ~b1) | ~b0 but it should be (~b2 & (~b1 | ~b0)) => (~b2 & ~b1) | (~b2 & ~b0))
  100   4   ~b2
  101   5                      (~b2 | ~b1) & ~b0
  110   6        ~b2 | ~b1
  111   7                      (~b2 | ~b1) | ~b0

 0000   0
 0001   1                                                  ~b3 & ~b2 & ~b1 & ~b0
 0010   2                          ~b3 & ~b2 & ~b1
 0011   3                                                  ~b3 & ~b2 & (~b1 | ~b0)
 0100   4            ~b3 & ~b2
 0101   5                                                  ~b3 & (~(~b2 | ~b1) & ~b0
 0110   6                          (~b3 & ~b2) | (~b3 & ~b1)
 0111   7                                                  (~b2 | ~b1) | ~b0
 1000   8   ~b3 | (true & ~b3)
 1001   9                                                       ~b3 | ((~b3 | ~b2) & ~b1 & ~b0)
 1010  10                          ~b3 | ((~b3 | ~b2) & ~b1)
 1011  11
 1100  12            ~b3 | (true & ~b2)
 1101  13
 1110  14                          ~b3 | ~b2 | ~b1
 1111  15                                                  ~b3 | ~b2 | ~b1 | ~b0


        */

        /** Heart of the algorithm
         *
         * A recursive algorithm checks if descending to left and/or right is required.
         * Aborts if any of the ranges have no characters. This accelerates computation on large alphabets.
         */
        auto rec = [&](this auto&& self, word_of_bits const& l_b1, word_of_bits const& r_b1, word_of_bits const& l_b2, word_of_bits const& r_b2, int level, size_t pr1_S, size_t pr2_S, size_t pr1_E, size_t pr2_E, size_t symb) -> void {
            auto lsymb = symb | (1 << level);
            auto b1 = l_b1 | (r_b1 & ~bits_lb[level]);
            auto b2 = l_b2 | (r_b2 & ~bits_rb[level]);

            auto [pr1, pr2] = count_pr(b1, b2, lsymb);

            assert(This->prefix_rank(idx1, lsymb) == pr1);
            assert(This->prefix_rank(idx2, lsymb) == pr2);

            auto rank_l_lb = pr1-pr1_S;
            auto rank_l_rb = pr2-pr2_S;

            auto rank_r_lb = pr1_E-pr1;
            auto rank_r_rb = pr2_E-pr2;
            if (level == 0) {
                assert(This->prefix_rank(idx1, symb) == pr1_S);
                assert(This->prefix_rank(idx2, symb) == pr2_S);
                cb( symb, rank_l_lb, rank_l_rb, pr1_S, pr2_S);
                cb(lsymb, rank_r_lb, rank_r_rb, pr1, pr2);
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

    /* This is a special variant of the `all_ranks_and_prefix_ranks` function.
     *
     * Differences:
     *  - it looks up the values for two positions at the same time.
     *  - it uses a callback function
     *  - it does not report (prefix) ranks for symbols that have 0-values
     */
    template <size_t L>
    void all_ranks_dual_limit(size_t idx1, size_t idx2, auto const& cb) const requires (L != 3) {
        assert(idx1 <= totalLength);
        assert(idx2 <= totalLength);

        // takes an index and splits it into the indices required for partial-paired-blocks
        auto split = [&](size_t idx) -> std::tuple<size_t, size_t, size_t, int64_t, int64_t> {
            auto bitId = idx % (l1_bits_ct*2);
            auto l1Id = idx / l1_bits_ct;
            auto l0Id = idx / l0_bits_ct;
            assert(l1Id < bits.size());
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
        auto const& bits_lb = bits[l1Id_lb].bits;
        auto const& bits_rb = bits[l1Id_rb].bits;

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
            auto b1 = l_b1 | (r_b1 & ~bits_lb[level]);
            auto b2 = l_b2 | (r_b2 & ~bits_rb[level]);

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

    template <size_t L>
    void all_ranks_dual_limit(size_t idx1, size_t idx2, auto const& cb) const requires (L == 3) {
        assert(idx1 <= totalLength);
        assert(idx2 <= totalLength);

        // takes an index and splits it into the indices required for partial-paired-blocks
        auto split = [&](size_t idx) -> std::tuple<size_t, size_t, size_t, int64_t, int64_t> {
            auto bitId = idx % (l1_bits_ct*2);
            auto l1Id = idx / l1_bits_ct;
            auto l0Id = idx / l0_bits_ct;
            assert(l1Id < bits.size());
            assert(l1Id/2 < l1.size());
            assert(l0Id/2 < l0.size());
            int64_t right_l1 = (l1Id%2)*2-1;
            int64_t right_l0 = (l0Id%2)*2-1;
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
        auto const& bits_lb = bits[l1Id_lb].bits;
        auto const& bits_rb = bits[l1Id_rb].bits;

        // computes prefix rank for symbol (symb)
        auto count_pr = [&](auto const& b1, auto const& b2, size_t symb) -> std::tuple<size_t, size_t> {
            auto count_lb = skip_first_or_last_n_bits_and_count(b1, bitId_lb);
            auto count_rb = skip_first_or_last_n_bits_and_count(b2, bitId_rb);

            auto pr1 = l0_lb[symb] + right_l0_lb * l1_lb[symb] + right_l1_lb * count_lb;
            auto pr2 = l0_rb[symb] + right_l0_rb * l1_rb[symb] + right_l1_rb * count_rb;
            assert(pr1 <= idx1);
            assert(pr2 <= idx2);

            return {pr1, pr2};
        };

        //using word_of_bits = std::bitset<l1_bits_ct>;
        auto This = this; //!TODO !WORKAROUND otherwise crashes gcc
        (void)This;

        // unroll
            //auto symb = 0;
            auto lsymb = 2;
            auto _b1 = ~bits_lb[1];
            auto _b2 = ~bits_rb[1];

            auto [pr1, pr2] = count_pr(_b1, _b2, lsymb);

            assert(This->prefix_rank(idx1, lsymb) == pr1);
            assert(This->prefix_rank(idx2, lsymb) == pr2);

            auto rank_l_lb = pr1-0;
            auto rank_l_rb = pr2-0;

            auto rank_r_lb = idx1-pr1;
            auto rank_r_rb = idx2-pr2;

            auto countLeft  = rank_l_rb - rank_l_lb;
            auto countRight = rank_r_rb - rank_r_lb;

            // recursion left
            if (countLeft > 0) {
                auto pr1_S = 0;
                auto pr2_S = 0;
                auto pr1_E = pr1;
                auto pr2_E = pr2;

                auto symb = 0;
                auto lsymb = 1;
                auto b1 = _b1 & ~bits_lb[0];
                auto b2 = _b2 & ~bits_rb[0];

                auto [pr1, pr2] = count_pr(b1, b2, lsymb);

                assert(This->prefix_rank(idx1, lsymb) == pr1);
                assert(This->prefix_rank(idx2, lsymb) == pr2);

                auto rank_l_lb = pr1-pr1_S;
                auto rank_l_rb = pr2-pr2_S;

                auto rank_r_lb = pr1_E-pr1;
                auto rank_r_rb = pr2_E-pr2;

                auto countLeft  = rank_l_rb - rank_l_lb;
                auto countRight = rank_r_rb - rank_r_lb;
                if (countLeft > 0) {
                    cb( symb >> (bitct-L), rank_l_lb, rank_l_rb, pr1_S, pr2_S);
                }
                if (countRight > 0) {
                    cb(lsymb >> (bitct-L), rank_r_lb, rank_r_rb, pr1, pr2);
                }
            }
            // recursion right
            if (countRight > 0) {
                auto pr1_S = pr1;
                auto pr2_S = pr2;
                auto pr1_E = idx1;
                auto pr2_E = idx2;

                auto symb = 2;
                auto lsymb = 3;
                auto b1 = _b1 | ~bits_lb[0];
                auto b2 = _b2 | ~bits_rb[0];

                auto [pr1, pr2] = count_pr(b1, b2, lsymb);

                assert(This->prefix_rank(idx1, lsymb) == pr1);
                assert(This->prefix_rank(idx2, lsymb) == pr2);

                auto rank_l_lb = pr1-pr1_S;
                auto rank_l_rb = pr2-pr2_S;

                auto rank_r_lb = pr1_E-pr1;
                auto rank_r_rb = pr2_E-pr2;

                auto countLeft  = rank_l_rb - rank_l_lb;
                auto countRight = rank_r_rb - rank_r_lb;
                if (countLeft > 0) {
                    cb( symb >> (bitct-L), rank_l_lb, rank_l_rb, pr1_S, pr2_S);
                }
                if (countRight > 0) {
                    cb(lsymb >> (bitct-L), rank_r_lb, rank_r_rb, pr1, pr2);
                }
            }
        // unroll end
    }



    template <typename Archive>
    void serialize(this auto&& self, Archive& ar) {
        ar(self.l0, self.l1, self.bits, self.totalLength);
    }
};

template <size_t Sigma> using PairedFlattenedBitvectors_64_4k   = PairedFlattenedBitvectors2L<Sigma, 64, 4096>;
template <size_t Sigma> using PairedFlattenedBitvectors_128_4k  = PairedFlattenedBitvectors2L<Sigma, 128, 4096>;
template <size_t Sigma> using PairedFlattenedBitvectors_256_4k  = PairedFlattenedBitvectors2L<Sigma, 256, 4096>;
template <size_t Sigma> using PairedFlattenedBitvectors_512_4k  = PairedFlattenedBitvectors2L<Sigma, 512, 4096>;
template <size_t Sigma> using PairedFlattenedBitvectors_1024_4k = PairedFlattenedBitvectors2L<Sigma, 1024, 4096>;
template <size_t Sigma> using PairedFlattenedBitvectors_2048_4k = PairedFlattenedBitvectors2L<Sigma, 2048, 4096>;

static_assert(checkString_c<PairedFlattenedBitvectors_64_4k>);
static_assert(checkString_c<PairedFlattenedBitvectors_128_4k>);
static_assert(checkString_c<PairedFlattenedBitvectors_256_4k>);
static_assert(checkString_c<PairedFlattenedBitvectors_512_4k>);
static_assert(checkString_c<PairedFlattenedBitvectors_1024_4k>);
static_assert(checkString_c<PairedFlattenedBitvectors_2048_4k>);

template <size_t Sigma> using PairedFlattenedBitvectors_64_64k   = PairedFlattenedBitvectors2L<Sigma, 64, 65536>;
template <size_t Sigma> using PairedFlattenedBitvectors_128_64k  = PairedFlattenedBitvectors2L<Sigma, 128, 65536>;
template <size_t Sigma> using PairedFlattenedBitvectors_256_64k  = PairedFlattenedBitvectors2L<Sigma, 256, 65536>;
template <size_t Sigma> using PairedFlattenedBitvectors_512_64k  = PairedFlattenedBitvectors2L<Sigma, 512, 65536>;
template <size_t Sigma> using PairedFlattenedBitvectors_1024_64k = PairedFlattenedBitvectors2L<Sigma, 1024, 65536>;
template <size_t Sigma> using PairedFlattenedBitvectors_2048_64k = PairedFlattenedBitvectors2L<Sigma, 2048, 65536>;
template <size_t Sigma> using PairedFlattenedBitvectors_4096_64k = PairedFlattenedBitvectors2L<Sigma, 4096, 65536>;

static_assert(checkString_c<PairedFlattenedBitvectors_64_64k>);
static_assert(checkString_c<PairedFlattenedBitvectors_128_64k>);
static_assert(checkString_c<PairedFlattenedBitvectors_256_64k>);
static_assert(checkString_c<PairedFlattenedBitvectors_512_64k>);
static_assert(checkString_c<PairedFlattenedBitvectors_1024_64k>);
static_assert(checkString_c<PairedFlattenedBitvectors_2048_64k>);
static_assert(checkString_c<PairedFlattenedBitvectors_4096_64k>);

template <size_t Sigma> using PairedFlattenedBitvectors_64_64kUA   = PairedFlattenedBitvectors2L<Sigma, 64, 65536, false>;
template <size_t Sigma> using PairedFlattenedBitvectors_128_64kUA  = PairedFlattenedBitvectors2L<Sigma, 128, 65536, false>;
template <size_t Sigma> using PairedFlattenedBitvectors_256_64kUA  = PairedFlattenedBitvectors2L<Sigma, 256, 65536, false>;
template <size_t Sigma> using PairedFlattenedBitvectors_512_64kUA  = PairedFlattenedBitvectors2L<Sigma, 512, 65536, false>;
template <size_t Sigma> using PairedFlattenedBitvectors_1024_64kUA = PairedFlattenedBitvectors2L<Sigma, 1024, 65536, false>;
template <size_t Sigma> using PairedFlattenedBitvectors_2048_64kUA = PairedFlattenedBitvectors2L<Sigma, 2048, 65536, false>;
template <size_t Sigma> using PairedFlattenedBitvectors_4096_64kUA = PairedFlattenedBitvectors2L<Sigma, 4096, 65536, false>;

static_assert(checkString_c<PairedFlattenedBitvectors_64_64kUA>);
static_assert(checkString_c<PairedFlattenedBitvectors_128_64kUA>);
static_assert(checkString_c<PairedFlattenedBitvectors_256_64kUA>);
static_assert(checkString_c<PairedFlattenedBitvectors_512_64kUA>);
static_assert(checkString_c<PairedFlattenedBitvectors_1024_64kUA>);
static_assert(checkString_c<PairedFlattenedBitvectors_2048_64kUA>);
static_assert(checkString_c<PairedFlattenedBitvectors_4096_64kUA>);

}
