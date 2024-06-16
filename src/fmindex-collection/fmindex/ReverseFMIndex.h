// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../occtable/concepts.h"
#include "../suffixarray/CSA.h"
#include "../utils.h"

namespace fmindex_collection {

/*
 * Same as the FMIndex, but build internally over the reverse text.
 *
 * This allows to have "extend_right" functionality instead of "extend_left".
 */
template <OccTable Table, SuffixArray_c TCSA = CSA>
struct ReverseFMIndex {
    static size_t constexpr Sigma = Table::Sigma;

    Table  occ;
    TCSA   csa;

    ReverseFMIndex() = default;
    ReverseFMIndex(std::span<uint8_t const> bwt, TCSA _csa)
        : occ{bwt}
        , csa{std::move(_csa)}
    {}

    ReverseFMIndex(Sequences auto const& _input, size_t samplingRate, size_t threadNbr) {

        auto [totalSize, inputText, inputSizes] = createSequences(_input, /*reverse*/ true);

        auto [bwt, csa] = [&, &inputText=inputText, &inputSizes=inputSizes] () {
            auto sa  = createSA(inputText, threadNbr);
            auto bwt = createBWT(inputText, sa);
            auto csa = TCSA{std::move(sa), samplingRate, inputSizes, /*reverse*/ true};

            return std::make_tuple(std::move(bwt), std::move(csa));
        }();

        decltype(inputText){}.swap(inputText); // inputText memory can be deleted

        *this = ReverseFMIndex{bwt, std::move(csa)};
    }

    size_t memoryUsage() const requires OccTableMemoryUsage<Table> {
        return occ.memoryUsage() + csa.memoryUsage();
    }

    size_t size() const {
        return occ.size();
    }

    auto locate(size_t idx) const -> std::tuple<size_t, size_t> {
        if constexpr (requires(Table t) {{ t.hasValue(size_t{}) }; }) {
            bool v = occ.hasValue(idx);
            uint64_t steps{};
            while(!v) {
                idx = occ.rank_symbol(idx);
                steps += 1;
                v = occ.hasValue(idx);
            }
            auto [chr, pos] = csa.value(idx);
            return {chr, pos-steps};

        } else {
            auto opt = csa.value(idx);
            uint64_t steps{};
            while(!opt) {
                if constexpr (requires(Table t) { { t.rank_symbol(size_t{}) }; }) {
                    idx = occ.rank_symbol(idx);
                } else {
                    idx = occ.rank(idx, occ.symbol(idx));
                }
                steps += 1;
                opt = csa.value(idx);
            }
            auto [chr, pos] = *opt;
            return {chr, pos-steps};
        }
    }

    auto single_locate_step(size_t idx) const -> std::optional<std::tuple<size_t, size_t>> {
        return csa.value(idx);
    }


    template <typename Archive>
    void serialize(Archive& ar) {
        ar(occ, csa);
    }
};

}
