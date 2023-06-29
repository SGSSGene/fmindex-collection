// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include "CSA_32.h"
#include "occtable/concepts.h"
#include "utils.h"

#include <algorithm>

namespace fmindex_collection {

template <OccTable_32 Table, typename TCSA = CSA_32>
struct BiFMIndex_32 {
    static size_t constexpr Sigma = Table::Sigma;

    using TTable = Table;

    Table  occ;
    Table  occRev;
    TCSA   csa;

//private:
    BiFMIndex_32(std::span<uint32_t const> bwt, std::vector<uint32_t> const& bwtRev, TCSA _csa)
        : occ{bwt}
        , occRev{bwtRev}
        , csa{std::move(_csa)}
    {
        assert(bwt.size() == bwtRev.size());
        assert(occ.size() == occRev.size());
        if (bwt.size() != bwtRev.size()) {
            throw std::runtime_error("bwt don't have the same size: " + std::to_string(bwt.size()) + " " + std::to_string(bwtRev.size()));
        }
        if (occ.size() != occRev.size()) {
            throw std::runtime_error("occ don't have the same size: " + std::to_string(occ.size()) + " " + std::to_string(occRev.size()));
        }
        // compute last row
        auto ct = std::array<uint64_t, Sigma>{};
        for (auto v : bwt) {
            ct[v] += 1;
        }
        for (size_t i{1}; i < ct.size(); ++i) {
            ct[i] = ct[i-1] + ct[i];
        }
        // check last row is correct
        for (size_t sym{0}; sym < Sigma; ++sym) {
            if (occ.rank(occ.size(), sym) != ct[sym]) {
                auto e = std::string{"Wrong rank for the last entry."}
                    + " Got different values for forward index."
                    + " sym: " + std::to_string(sym)
                    + " got: " + std::to_string(occ.rank(occ.size(), sym))
                    + " expected: " + std::to_string(ct[sym]);
                throw std::runtime_error(e);
            }
            if (occRev.rank(occRev.size(), sym) != ct[sym]) {
                auto e = std::string{"Wrong rank for the last entry."}
                    + " Got different values for reverse index."
                    + " sym: " + std::to_string(sym)
                    + " got: " + std::to_string(occRev.rank(occRev.size(), sym))
                    + " expected: " + std::to_string(ct[sym]);
                throw std::runtime_error(e);
            }
        }
        if constexpr (requires(Table t) {{ t.hasValue(size_t{}) }; }) {
            for (size_t i{0}; i < occ.size(); ++i) {
                if (csa.value(i).has_value()) {
                    occ.setValue(i);
                }
            }
        }
    }

public:
    /**!\brief Creates a BiFMIndex with a specified sampling rate
     *
     * \param _input a list of sequences
     * \param samplingRate rate of the sampling
     */
    BiFMIndex_32(Sequences_32 auto const& _input, size_t samplingRate, size_t threadNbr)
        : occ{cereal_tag{}}
        , occRev{cereal_tag{}}
        , csa{cereal_tag{}}
    {
        auto [totalSize, inputText, inputSizes] = createSequences_32(_input, samplingRate);

        // create BurrowsWheelerTransform and CompressedSuffixArray
        auto [bwt, csa] = [&, &inputText=inputText, &inputSizes=inputSizes] () {
            auto sa  = createSA_32(inputText, threadNbr);
            auto bwt = createBWT_32(inputText, sa);
            auto csa = TCSA(std::move(sa), samplingRate, inputSizes);
            return std::make_tuple(std::move(bwt), std::move(csa));
        }();

        // create BurrowsWheelerTransform on reversed text
        auto bwtRev = [&, &inputText=inputText]() {
            std::ranges::reverse(inputText);
            auto saRev  = createSA_32(inputText, threadNbr);
            auto bwtRev = createBWT_32(inputText, saRev);
            return bwtRev;
        }();

        decltype(inputText){}.swap(inputText); // inputText memory can be deleted

        *this = BiFMIndex_32{bwt, bwtRev, std::move(csa)};
    }


    /*!\brief Specific c'tor for serialization use
     */
    BiFMIndex_32(cereal_tag)
        : occ{cereal_tag{}}
        , occRev{cereal_tag{}}
        , csa{cereal_tag{}}
    {}

    size_t memoryUsage() const requires OccTableMemoryUsage<Table> {
        return occ.memoryUsage() + occRev.memoryUsage() + csa.memoryUsage();
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
            return {chr, pos+steps};

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
            return {chr, pos+steps};
        }
    }

    auto locate(size_t idx, size_t maxSteps) const -> std::optional<std::tuple<size_t, size_t>> {
        auto opt = csa.value(idx);
        uint64_t steps{};
        for (;!opt and maxSteps > 0; --maxSteps) {
            idx = occ.rank(idx, occ.symbol(idx));
            steps += 1;
            opt = csa.value(idx);
        }
        if (opt) {
            std::get<1>(*opt) += steps;
        }
        return opt;
    }


    auto single_locate_step(size_t idx) const -> std::optional<std::tuple<size_t, size_t>> {
        return csa.value(idx);
    }


    template <typename Archive>
    void serialize(Archive& ar) {
        ar(occ, occRev, csa);
    }
};

}
