// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../string/concepts.h"
#include "../suffixarray/CSA.h"
#include "../utils.h"

namespace fmindex_collection {

/*
 * Same as the FMIndex, but build internally over the reverse text.
 *
 * This allows to have "extend_right" functionality instead of "extend_left".
 */
template <String_c String, SuffixArray_c TCSA = CSA>
struct ReverseFMIndex {
    using ADEntry = std::tuple<size_t, size_t>;
    static size_t constexpr Sigma = String::Sigma;

    String                      bwt;
    std::array<size_t, Sigma+1> C{};
    TCSA   csa;

    ReverseFMIndex() = default;
    ReverseFMIndex(std::span<uint8_t const> _bwt, TCSA _csa)
        : bwt{_bwt}
        , csa{std::move(_csa)}
    {
        for (auto c : _bwt) {
            C[c+1] += 1;
        }
        for (size_t i{1}; i < C.size(); ++i) {
            C[i] = C[i] + C[i-1];
        }
    }

    ReverseFMIndex(Sequences auto const& _input, size_t samplingRate, size_t threadNbr) {

        auto [totalSize, inputText, inputSizes] = createSequences(_input, /*reverse*/ true);

        auto [bwt, csa] = [&, &inputText=inputText, &inputSizes=inputSizes] () {
            auto sa  = createSA64(inputText, threadNbr);
            auto bwt = createBWT64(inputText, sa);
            auto csa = TCSA{std::move(sa), samplingRate, inputSizes, /*reverse*/ true};

            return std::make_tuple(std::move(bwt), std::move(csa));
        }();

        decltype(inputText){}.swap(inputText); // inputText memory can be deleted

        *this = ReverseFMIndex{bwt, std::move(csa)};
    }

/*    size_t memoryUsage() const requires OccTableMemoryUsage<Table> {
        return occ.memoryUsage() + csa.memoryUsage();
    }*/

    size_t size() const {
        return bwt.size();
    }

    auto locate(size_t idx) const -> std::tuple<ADEntry, size_t> {
        if constexpr (requires(String t) {{ t.hasValue(size_t{}) }; }) {
            bool v = bwt.hasValue(idx);
            uint64_t steps{};
            while(!v) {
                idx = bwt.rank_symbol(idx);
                steps += 1;
                v = bwt.hasValue(idx);
            }
            //!TODO steps is from the end???
            return {csa.value(idx), steps};

        } else {
            auto opt = csa.value(idx);
            uint64_t steps{};
            while(!opt) {
                if constexpr (requires(String t) { { t.rank_symbol(size_t{}) }; }) {
                    idx = bwt.rank_symbol(idx);
                } else {
                    auto symb = bwt.symbol(idx);
                    idx = bwt.rank(idx, symb) + C[symb];
                }
                steps += 1;
                opt = csa.value(idx);
            }
            //!TODO steps is from the end???
            return {*opt, steps};
        }
    }

    auto single_locate_step(size_t idx) const -> std::optional<ADEntry> {
        return csa.value(idx);
    }


    template <typename Archive>
    void serialize(Archive& ar) {
        ar(bwt, C, csa);
    }
};

}
