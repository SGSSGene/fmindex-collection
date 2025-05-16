// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../string/concepts.h"
#include "../string/utils.h"
#include "../suffixarray/CSA.h"
#include "../utils.h"

namespace fmindex_collection {

template <String_c String, SuffixArray_c TCSA = CSA>
struct FMIndex {
    static size_t constexpr Sigma = String::Sigma;

    String                      bwt;
    std::array<size_t, Sigma+1> C{0};
    TCSA   csa;

    FMIndex() = default;
    FMIndex(FMIndex const&) = delete;
    FMIndex(FMIndex&&) noexcept = default;
    FMIndex(std::span<uint8_t const> _bwt, TCSA _csa)
        : bwt{_bwt}
        , C{computeAccumulatedC(bwt)}
        , csa{std::move(_csa)}
    {}

    FMIndex(std::vector<uint8_t> _input, size_t samplingRate, size_t threadNbr) {
        auto input = std::vector<std::vector<uint8_t>>{std::move(_input)};
        *this = FMIndex{std::move(input), samplingRate, threadNbr};
    }

    FMIndex(Sequences auto const& _input, size_t samplingRate, size_t threadNbr) {
        auto [totalSize, inputText, inputSizes] = createSequences(_input);

        if (totalSize < std::numeric_limits<int32_t>::max()) { // only 32bit SA required
            auto [bwt, csa] = [&]() {
                auto sa  = createSA32(inputText, threadNbr);
                auto bwt = createBWT32(inputText, sa);
                auto csa = TCSA{std::move(sa), samplingRate, inputSizes};

                return std::make_tuple(std::move(bwt), std::move(csa));
            }();

            *this = FMIndex{bwt, std::move(csa)};

        } else { // required 64bit SA required
            auto [bwt, csa] = [&]() {
                auto sa  = createSA64(inputText, threadNbr);
                auto bwt = createBWT64(inputText, sa);
                auto csa = TCSA{std::move(sa), samplingRate, inputSizes};

                return std::make_tuple(std::move(bwt), std::move(csa));
            }();

            *this = FMIndex{bwt, std::move(csa)};
        }
    }
    auto operator=(FMIndex const&) -> FMIndex& = delete;
    auto operator=(FMIndex&&) noexcept -> FMIndex& = default;


/*    size_t memoryUsage() const requires OccTableMemoryUsage<Table> {
        return occ.memoryUsage() + csa.memoryUsage();
    }*/

    size_t size() const {
        return bwt.size();
    }

    auto locate(size_t idx) const -> std::tuple<size_t, size_t> {
        auto opt = csa.value(idx);
        size_t steps{};
        while(!opt) {
            auto symb = bwt.symbol(idx);
            idx = bwt.rank(idx, symb) + C[symb];
            steps += 1;
            opt = csa.value(idx);
        }
        auto [chr, pos] = *opt;
        return std::make_tuple(chr, pos + steps);
    }

    auto single_locate_step(size_t idx) const -> std::optional<std::tuple<size_t, size_t>> {
        return csa.value(idx);
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(bwt, C, csa);
    }
};

}
