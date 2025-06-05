// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../builtins.h"
#include "../bitvector/Bitvector.h"
#include "../string/PairedFlattenBitvectors_L0L1.h"
#include "concepts.h"

#include <algorithm>
#include <array>
#include <bit>
#include <bitset>
#include <cassert>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <vector>

namespace fmindex_collection::string {

/* Implements the concept `String_c`
 *
 * \param TSigma size of the alphabet
 */
template <size_t TSigma, template <size_t> typename String_t0 = PairedFlattenBitvectors_512_64k, size_t L0Size = 1ull<<(std::bit_width(TSigma)/2),
                         template <size_t> typename String_t1 = String_t0, size_t L1Size = std::max((TSigma+L0Size-1)/L0Size, size_t{2})>
struct MultiaryWavelet {
    using StringL0 = String_t0<L0Size>;
    using StringL1 = String_t1<L1Size>;

    static constexpr size_t Sigma = TSigma;

    StringL0                     l0;
    std::array<StringL1, L0Size> l1;

    static_assert(TSigma <= StringL0::Sigma * StringL1::Sigma);


    MultiaryWavelet()
        : MultiaryWavelet{internal_tag{}, std::span<uint8_t const>{}}
    {}

    template <typename = uint8_t>
    MultiaryWavelet(std::span<uint8_t const> _symbols)
        : MultiaryWavelet{internal_tag{}, _symbols}
    {
        static_assert(Sigma < 256, "This constructor can only be used, if Alphabet size is smaller than 256");
    }

    MultiaryWavelet(std::span<uint64_t const> _symbols)
        : MultiaryWavelet{internal_tag{}, _symbols}
    {}

private:
    struct internal_tag{};
    template <typename T>
    MultiaryWavelet(internal_tag, std::span<T const> _symbols) {
        //!TODO can we do this without std::vector?
        auto l0Buffer = std::vector<T>{};
        auto l1Buffers = std::array<std::vector<T>, StringL0::Sigma>{};
        for (auto c : _symbols) {
            auto l0c = c / StringL1::Sigma;
            l0Buffer.push_back(l0c);
            auto l1c = c % StringL1::Sigma;
            l1Buffers[l0c].push_back(l1c);
        }

        l0 = StringL0{l0Buffer};
        for (size_t i{0}; i < StringL0::Sigma; ++i) {
            l1[i] = StringL1{l1Buffers[i]};
        }
    }
public:

    size_t size() const {
        return l0.size();
    }

    uint64_t symbol(uint64_t idx) const {
        assert(idx < size());
        auto l0c = l0.symbol(idx);
        auto r   = l0.rank(idx, l0c);

        auto l1c = l1[l0c].symbol(r);
        return l0c * StringL1::Sigma + l1c;
    }

    uint64_t rank(uint64_t idx, uint64_t symb) const {
        assert(idx <= size());
        assert(symb < TSigma);

        auto l0c = symb / StringL1::Sigma;
        auto l1c = symb % StringL1::Sigma;
        auto r = l0.rank(idx, l0c);
        return l1[l0c].rank(r, l1c);
    }

    uint64_t prefix_rank(uint64_t idx, uint64_t symb) const {
        assert(idx <= size());
        assert(symb <= TSigma);

        auto l0c = symb / StringL1::Sigma;
        auto l1c = symb % StringL1::Sigma;
        auto pr = l0.prefix_rank(idx, l0c);
        auto r  = l0.rank(idx, l0c);
        return pr + l1[l0c].prefix_rank(r, l1c);
    }

    auto all_ranks(uint64_t idx) const -> std::array<uint64_t, TSigma> {
        assert(idx <= size());
        auto r = std::array<uint64_t, TSigma>{};
        for (size_t symb{0}; symb < TSigma; ++symb) {
            r[symb] = rank(idx, symb);
        }
        return r;
    }

    auto all_ranks_and_prefix_ranks(uint64_t idx) const -> std::tuple<std::array<uint64_t, TSigma>, std::array<uint64_t, TSigma>> {
        assert(idx <= size());
        auto rs = all_ranks(idx);
        auto prs = std::array<uint64_t, TSigma>{};
        for (size_t i{1}; i < prs.size(); ++i) {
            prs[i] = prs[i-1] + rs[i-1];
        }
        return {rs, prs};
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(l0, l1);
    }
};

template <size_t Sigma> using MultiaryWavelet_64_64k  = MultiaryWavelet<Sigma, PairedFlattenBitvectors_64_64k>;
template <size_t Sigma> using MultiaryWavelet_512_64k = MultiaryWavelet<Sigma, PairedFlattenBitvectors_512_64k>;
template <size_t Sigma> using MultiaryWavelet_s16     = MultiaryWavelet<Sigma, PairedFlattenBitvectors_512_64k, 4, MultiaryWavelet>;
template <size_t Sigma> using MultiaryWavelet_s256    = MultiaryWavelet<Sigma, PairedFlattenBitvectors_512_64k, 8, MultiaryWavelet>;

static_assert(checkString_c<MultiaryWavelet_64_64k>);
static_assert(checkString_c<MultiaryWavelet_512_64k>);
static_assert(checkString_c<MultiaryWavelet_s16>);
static_assert(checkString_c<MultiaryWavelet_s256>);

}
