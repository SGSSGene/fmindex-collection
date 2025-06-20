// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "concepts.h"

#include "../bitvector/concepts.h"
#include "../bitvector/L0L1_NBitvector.h"

#include <bitset>
#include <cassert>
#include <ranges>
#include <vector>

namespace fmc::string {

template <Bitvector_c Bitvector>
struct WrappedBitvector {
    static constexpr size_t Sigma = 2;

    Bitvector bitvector;

    WrappedBitvector(std::span<uint8_t const> _symbols)
        : bitvector{ _symbols | std::views::transform([&](size_t s) -> bool {
            return s == 1;
        })}
    {}

    WrappedBitvector() = default;
    WrappedBitvector(WrappedBitvector const&) = default;
    WrappedBitvector(WrappedBitvector&&) noexcept = default;
    auto operator=(WrappedBitvector const&) -> WrappedBitvector& = default;
    auto operator=(WrappedBitvector&&) noexcept -> WrappedBitvector& = default;


    size_t size() const noexcept {
        return bitvector.size();
    }

    uint8_t symbol(size_t idx) const noexcept {
        return bitvector.symbol(idx);
    }

    uint64_t rank(size_t idx, uint8_t symb) const noexcept {
        auto r = bitvector.rank(idx);
//        return (idx - r) * (1-symb) + r * symb;
        return (idx - r + symb*(2*r-idx));

//        idx


/*        idx-r -symb(idx-r)    + r*symb;
        idx-r -symb*idx + symb*r + r*symb;

        idx-r-symb*idx + 2*symb*r

        idx(-r - symb) + 2*symb*r*/

/*        idx - r + r*symb + r*symb;
        idx - r + 2 * r * symb
        idx + 2*r*symb - r
        idx + r * (2*symb - 1);*/




//        return idx * (1-symb) - ((symb*2)-1) * r;

//        idx - idx*symb - ((symb*2*r) - r)
//        idx - idx*symb - 2*symb*r + r
//        idx + r - idx*symb - 2*symb*r
//        idx + r - symb * (idx + 2*r)

/*        if (symb) return r;
        return idx - r;*/
    }

    uint64_t prefix_rank(size_t idx, uint8_t symb) const noexcept {
        if (symb == 0) return 0;
        if (symb == 2) return idx;
        return idx - bitvector.rank(idx);
    }

    auto all_ranks(size_t idx) const -> std::array<uint64_t, Sigma> {
        auto r = bitvector.rank(idx);
        return {idx-r, r};
    }

    auto all_ranks_and_prefix_ranks(size_t idx) const -> std::tuple<std::array<uint64_t, Sigma>, std::array<uint64_t, Sigma>> {
        assert(idx <= size());

        auto rs  = all_ranks(idx);
        auto prs = std::array<uint64_t, Sigma>{};
        for (size_t i{1}; i < prs.size(); ++i) {
            prs[i] = prs[i-1] + rs[i-1];
        }
        return {rs, prs};
    }


    template <typename Archive>
    void serialize(Archive& ar) {
        ar(bitvector);
    }
};

static_assert(String_c<WrappedBitvector<bitvector::L0L1_512_64kBitvector>>);

}
