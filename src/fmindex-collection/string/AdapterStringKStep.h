// SPDX-FileCopyrightText: 2026 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "concepts.h"

#include <bit>
#include <ranges>

/** Takes a normal string that fulfills String_c concept
 * and extends it such that it fulfills the StringKStep_c concept
 */
namespace fmc::string {

template <size_t TSigma, size_t TL, template <size_t> typename TString>
    requires String_c<TString<TSigma>>
        && String_c<TString<(size_t{1}<<TL)>>
struct AdapterStringKStep {

    static constexpr size_t Sigma = TSigma;
    static constexpr size_t LSigma = (size_t{1}<<TL);
    static constexpr size_t L = TL;

    // number of full length bit vectors needed `2^bitct > TSigma`
    static constexpr auto bitct = std::bit_width(Sigma-1);

    TString<Sigma>  string;
    TString<LSigma> string_limit;

    AdapterStringKStep() = default;
    AdapterStringKStep(AdapterStringKStep&&) = default;
    auto operator=(AdapterStringKStep&&) -> AdapterStringKStep& = default;


    template <std::ranges::range range_t>
        requires std::convertible_to<std::ranges::range_value_t<range_t>, uint64_t>
    AdapterStringKStep(range_t&& _symbols)
        : string{_symbols}
        , string_limit{ _symbols | std::views::transform([](uint64_t v) {
            return (v >> (bitct - TL));
        })}
    {
        assert(string.size() == string_limit.size());
    }

    auto size() const {
        return string.size();
    }
    auto symbol(size_t idx) const {
        return string.symbol(idx);
    }
    auto rank(size_t idx, uint64_t symb) const {
        return string.rank(idx, symb);
    }
    auto prefix_rank(size_t idx, uint64_t symb) const {
        return string.prefix_rank(idx, symb);
    }
    auto all_ranks(size_t idx) const {
        return string.all_ranks(idx);
    }
    auto all_ranks_and_prefix_ranks(size_t idx) const {
        return string.all_ranks_and_prefix_ranks(idx);
    }

    template <size_t L2 = TL>
    auto symbol_limit(size_t idx) const {
        if constexpr (L2 == bitct) {
            return string.symbol(idx);
        } else if constexpr (L2 == TL) {
            return string_limit.symbol(idx);
        } else {
            []<bool b = false>() {
                static_assert(b, "Only works for L2 is either TL or bitct");
            }();
        }
    }

    template <size_t L2 = TL>
    auto rank_limit(size_t idx, uint64_t symb) const {
        if constexpr (L2 == bitct) {
            return string.rank(idx, symb);
        } else if constexpr (L2 == TL) {
            return string_limit.rank(idx, symb);
        } else {
            []<bool b = false>() {
                static_assert(b, "Only works for L2 is either TL or bitct");
            }();
        }
    }

    template <size_t L2 = TL>
    auto prefix_rank_and_rank_limit(size_t idx, uint64_t symb) const -> std::tuple<uint64_t, uint64_t> {
        if constexpr (L2 == bitct) {
            auto r  = string.rank(idx, symb);
            auto pr = string.prefix_rank(idx, symb);
            return {pr, r};
        } else if constexpr (L2 == TL) {
            auto r  = string_limit.rank(idx, symb);
            auto pr = string_limit.prefix_rank(idx, symb);
            return {pr, r};
        } else {
            []<bool b = false>() {
                static_assert(b, "Only works for L2 is either TL or bitct");
            }();
        }
    }


    template <typename Archive>
    void serialize(this auto&& self, Archive& ar) {
        ar(self.string, self.string_limit);
    }
};

}
