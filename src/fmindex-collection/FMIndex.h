// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include "CSA.h"
#include "occtable/concepts.h"
#include "utils.h"

namespace fmindex_collection {

template <OccTable Table, typename TCSA = CSA>
struct FMIndex {
    static size_t constexpr Sigma = Table::Sigma;

    using TTable = Table;

    Table  occ;
    TCSA   csa;

    FMIndex(std::span<uint8_t const> bwt, TCSA _csa)
        : occ{bwt}
        , csa{std::move(_csa)}
    {}

    FMIndex(std::vector<uint8_t> _input, size_t samplingRate, size_t threadNbr)
        : occ{cereal_tag{}}
        , csa{cereal_tag{}}
    {
        auto input = std::vector<std::vector<uint8_t>>{std::move(_input)};
        *this = FMIndex{std::move(input), samplingRate, threadNbr};
    }

    FMIndex(Sequences auto const& _input, size_t samplingRate, size_t threadNbr)
        : occ{cereal_tag{}}
        , csa{cereal_tag{}}
    {
        auto [totalSize, inputText, inputSizes] = createSequences(_input, samplingRate);

        auto [bwt, csa] = [&, &inputText=inputText, &inputSizes=inputSizes] () {
            auto sa  = createSA(inputText, threadNbr);
            auto bwt = createBWT(inputText, sa);
            auto csa = TCSA{std::move(sa), samplingRate, inputSizes};

            return std::make_tuple(std::move(bwt), std::move(csa));
        }();

        *this = FMIndex{bwt, std::move(csa)};
    }


    FMIndex(cereal_tag)
        : occ{cereal_tag{}}
        , csa{cereal_tag{}}
    {}

    size_t memoryUsage() const requires OccTableMemoryUsage<Table> {
        return occ.memoryUsage() + csa.memoryUsage();
    }

    size_t size() const {
        return occ.size();
    }

    auto locate(size_t idx) const -> std::tuple<size_t, size_t> {
        auto opt = csa.value(idx);
        uint64_t steps{};
        while(!opt) {
            idx = occ.rank(idx, occ.symbol(idx));
            steps += 1;
            opt = csa.value(idx);
        }
        auto [chr, pos] = *opt;
        return std::make_tuple(chr, pos + steps);
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
        ar(occ, csa);
    }
};

}
