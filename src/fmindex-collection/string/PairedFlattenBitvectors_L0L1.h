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


namespace fmindex_collection::string {


template <size_t TSigma, size_t l1_bits_ct, size_t l0_bits_ct, bool Align=true>
struct PairedFlattenBitvectors_L0L1 {
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
            auto v = neprv8_detail::rank(bits, symb);
            return skip_first_or_last_n_bits_and_count(v, idx);
        }

        uint64_t prefix_rank(uint64_t idx, uint64_t symb) const {
            assert(idx <= l1_bits_ct*2);
            assert(symb <= Sigma);
            auto v = neprv8_detail::prefix_rank(bits, symb);
            return skip_first_or_last_n_bits_and_count(v, idx);
        }

        auto all_ranks(uint64_t idx) const -> std::array<uint64_t, TSigma> {
            assert(idx <= l1_bits_ct*2);

            auto vs = neprv8_detail::rank_all<(1ull<<bitct)>(bits);
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

    using BlockL1 = std::array<uint16_t, TSigma+1>;
    using BlockL0 = std::array<uint64_t, TSigma+1>;

    std::vector<InBits> bits{{}};
    std::vector<BlockL1> l1{{}};
    std::vector<BlockL0> l0{{}};
    size_t totalLength{};

    PairedFlattenBitvectors_L0L1() = default;

    PairedFlattenBitvectors_L0L1(std::span<uint8_t const> _symbols)
        : PairedFlattenBitvectors_L0L1{internal_tag{}, _symbols}
    {}

    PairedFlattenBitvectors_L0L1(std::span<uint64_t const> _symbols)
        : PairedFlattenBitvectors_L0L1{internal_tag{}, _symbols}
    {}

private:
    struct internal_tag{};
    template <typename T>
    PairedFlattenBitvectors_L0L1(internal_tag, std::span<T const> _symbols) {
        auto const _length = _symbols.size();
        bits.reserve(_length/l1_bits_ct + 2);
        if (_length == 0) return;

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
                        acc[symb+1] += a;
                    }
                    if (i % 2 == 0) {
                        for (size_t symb{0}; symb <= TSigma; ++symb) {
                            l1[l0I*l1_block_ct + i/2][symb] = acc[symb];
                        }
                    }
                }
                for (size_t symb{0}; symb <= TSigma; ++symb) {
                    l0_acc[symb] += acc[symb];
                }
                // update l0 (reached center)
                l0[l0I] = l0_acc;
                // walk backwards through left part and revert l0
                for (size_t i{0}; i < l1_block_ct; ++i) {
                    if (i % 2 == 0) {
                        for (size_t symb{0}; symb <= TSigma; ++symb) {
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
                        acc[symb+1] += a;
                    }
                    if (i % 2 == 0) {
                        for (size_t symb{0}; symb <= TSigma; ++symb) {
                            l1[l0I*l1_block_ct + i/2][symb] = acc[symb];
                        }
                    }
                }

                for (size_t symb{0}; symb <= TSigma; ++symb) {
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
        assert(l1Id < bits.size());
        auto symb = bits[l1Id].symbol(bitId);
        assert(symb < Sigma);
        return symb;
    }

    uint64_t rank(uint64_t idx, uint64_t symb) const {
        assert(idx <= totalLength);
        assert(symb < Sigma);
        auto bitId = idx % (l1_bits_ct*2);
        auto l1Id = idx / l1_bits_ct;
        auto l0Id = idx / l0_bits_ct;
        assert(l1Id < bits.size());
        assert(l1Id/2 < l1.size());
        assert(l0Id/2 < l0.size());

        int64_t right_l1 = (l1Id%2)*2-1;
        int64_t right_l0 = (l0Id%2)*2-1;

        auto count = bits[l1Id].rank(bitId, symb);

        auto r = (l0[l0Id/2][symb+1] - l0[l0Id/2][symb]) + right_l0 * (l1[l1Id/2][symb+1] - l1[l1Id/2][symb]) + right_l1 * count;
        assert(r <= idx);
        return r;
    }

    uint64_t prefix_rank(uint64_t idx, uint64_t symb) const {
        assert(idx <= totalLength);
        assert(symb <= Sigma);
        auto bitId = idx % (l1_bits_ct*2);
        auto l1Id = idx / l1_bits_ct;
        auto l0Id = idx / l0_bits_ct;
        assert(l1Id < bits.size());
        assert(l1Id/2 < l1.size());
        assert(l0Id/2 < l0.size());

        int64_t right_l1 = (l1Id%2)*2-1;
        int64_t right_l0 = (l0Id%2)*2-1;

        auto count = bits[l1Id].prefix_rank(bitId, symb);
        auto r = l0[l0Id/2][symb] + right_l0 * l1[l1Id/2][symb] + right_l1 * count;
        assert(r <= idx);
        return r;
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
    void serialize(Archive& ar) {
        ar(l0, l1, bits, totalLength);
    }
};

template <size_t Sigma> using PairedFlattenBitvectors_64_4k   = PairedFlattenBitvectors_L0L1<Sigma, 64, 4096>;
template <size_t Sigma> using PairedFlattenBitvectors_128_4k  = PairedFlattenBitvectors_L0L1<Sigma, 128, 4096>;
template <size_t Sigma> using PairedFlattenBitvectors_256_4k  = PairedFlattenBitvectors_L0L1<Sigma, 256, 4096>;
template <size_t Sigma> using PairedFlattenBitvectors_512_4k  = PairedFlattenBitvectors_L0L1<Sigma, 512, 4096>;
template <size_t Sigma> using PairedFlattenBitvectors_1024_4k = PairedFlattenBitvectors_L0L1<Sigma, 1024, 4096>;
template <size_t Sigma> using PairedFlattenBitvectors_2048_4k = PairedFlattenBitvectors_L0L1<Sigma, 2048, 4096>;

static_assert(checkString_c<PairedFlattenBitvectors_64_4k>);
static_assert(checkString_c<PairedFlattenBitvectors_128_4k>);
static_assert(checkString_c<PairedFlattenBitvectors_256_4k>);
static_assert(checkString_c<PairedFlattenBitvectors_512_4k>);
static_assert(checkString_c<PairedFlattenBitvectors_1024_4k>);
static_assert(checkString_c<PairedFlattenBitvectors_2048_4k>);

template <size_t Sigma> using PairedFlattenBitvectors_64_64k   = PairedFlattenBitvectors_L0L1<Sigma, 64, 65536>;
template <size_t Sigma> using PairedFlattenBitvectors_128_64k  = PairedFlattenBitvectors_L0L1<Sigma, 128, 65536>;
template <size_t Sigma> using PairedFlattenBitvectors_256_64k  = PairedFlattenBitvectors_L0L1<Sigma, 256, 65536>;
template <size_t Sigma> using PairedFlattenBitvectors_512_64k  = PairedFlattenBitvectors_L0L1<Sigma, 512, 65536>;
template <size_t Sigma> using PairedFlattenBitvectors_1024_64k = PairedFlattenBitvectors_L0L1<Sigma, 1024, 65536>;
template <size_t Sigma> using PairedFlattenBitvectors_2048_64k = PairedFlattenBitvectors_L0L1<Sigma, 2048, 65536>;
template <size_t Sigma> using PairedFlattenBitvectors_4096_64k = PairedFlattenBitvectors_L0L1<Sigma, 4096, 65536>;

static_assert(checkString_c<PairedFlattenBitvectors_64_64k>);
static_assert(checkString_c<PairedFlattenBitvectors_128_64k>);
static_assert(checkString_c<PairedFlattenBitvectors_256_64k>);
static_assert(checkString_c<PairedFlattenBitvectors_512_64k>);
static_assert(checkString_c<PairedFlattenBitvectors_1024_64k>);
static_assert(checkString_c<PairedFlattenBitvectors_2048_64k>);
static_assert(checkString_c<PairedFlattenBitvectors_4096_64k>);

template <size_t Sigma> using PairedFlattenBitvectors_64_64kUA   = PairedFlattenBitvectors_L0L1<Sigma, 64, 65536, false>;
template <size_t Sigma> using PairedFlattenBitvectors_128_64kUA  = PairedFlattenBitvectors_L0L1<Sigma, 128, 65536, false>;
template <size_t Sigma> using PairedFlattenBitvectors_256_64kUA  = PairedFlattenBitvectors_L0L1<Sigma, 256, 65536, false>;
template <size_t Sigma> using PairedFlattenBitvectors_512_64kUA  = PairedFlattenBitvectors_L0L1<Sigma, 512, 65536, false>;
template <size_t Sigma> using PairedFlattenBitvectors_1024_64kUA = PairedFlattenBitvectors_L0L1<Sigma, 1024, 65536, false>;
template <size_t Sigma> using PairedFlattenBitvectors_2048_64kUA = PairedFlattenBitvectors_L0L1<Sigma, 2048, 65536, false>;
template <size_t Sigma> using PairedFlattenBitvectors_4096_64kUA = PairedFlattenBitvectors_L0L1<Sigma, 4096, 65536, false>;

static_assert(checkString_c<PairedFlattenBitvectors_64_64kUA>);
static_assert(checkString_c<PairedFlattenBitvectors_128_64kUA>);
static_assert(checkString_c<PairedFlattenBitvectors_256_64kUA>);
static_assert(checkString_c<PairedFlattenBitvectors_512_64kUA>);
static_assert(checkString_c<PairedFlattenBitvectors_1024_64kUA>);
static_assert(checkString_c<PairedFlattenBitvectors_2048_64kUA>);
static_assert(checkString_c<PairedFlattenBitvectors_4096_64kUA>);

}
