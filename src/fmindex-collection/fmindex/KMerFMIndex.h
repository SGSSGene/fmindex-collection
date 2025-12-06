// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../bitvector/Bitvector2L.h"
#include "../string/concepts.h"
#include "../string/utils.h"
#include "../suffixarray/CSA.h"
#include "../utils.h"

namespace fmc {

template <String_c String, size_t KMer, SuffixArray_c TCSA = CSA>
struct KMerFMIndex {
    using ADEntry = std::tuple<size_t, size_t>;
    static size_t constexpr Sigma = String::Sigma;

    String                      bwt;
    std::array<size_t, Sigma+1> C{0};
    TCSA                        csa;
    bitvector::Bitvector2L<512, 65536> kmerStarts;


    KMerFMIndex() = default;
    KMerFMIndex(KMerFMIndex const&) = delete;
    KMerFMIndex(KMerFMIndex&&) noexcept = default;
    KMerFMIndex(std::span<uint8_t const> _bwt, TCSA _csa)
        : bwt{_bwt}
        , C{computeAccumulatedC(bwt)}
        , csa{std::move(_csa)}
    {
        // compute kmer start position, do a partial backwards search
        auto bits = std::vector<bool>{};
        bits.resize(bwt.size()+1);
        auto step = std::function<void(size_t, size_t, size_t)>{};
        step = [&](size_t lb, size_t rb, size_t depth) {
            bits[lb] = true;
            bits[rb] = true;
            if (lb == rb) return;
            if (depth >= KMer) return;
            for (size_t symb{0}; symb < Sigma; ++symb) {
                auto nlb = C[symb] + bwt.rank(lb, symb);
                auto nrb = C[symb] + bwt.rank(rb, symb);
                step(nlb, nrb, depth+1);
            }
        };
        step(0, bwt.size(), 0);
        kmerStarts = {bits};
    }

    KMerFMIndex(std::vector<uint8_t> _input, size_t samplingRate, size_t threadNbr, size_t seqOffset=0) {
        auto input = std::vector<std::vector<uint8_t>>{std::move(_input)};
        *this = KMerFMIndex{std::move(input), samplingRate, threadNbr, seqOffset};
    }

    KMerFMIndex(Sequences auto const& _input, size_t samplingRate, size_t threadNbr, size_t seqOffset=0) {
        auto [totalSize, inputText, inputSizes] = createSequences(_input);

        if (totalSize < std::numeric_limits<int32_t>::max()) { // only 32bit SA required
            auto [bwt, csa] = [&]() {
                auto sa  = createSA32(inputText, threadNbr);
                auto bwt = createBWT32(inputText, sa);
                auto csa = TCSA{std::move(sa), samplingRate, inputSizes, /*.reverse=*/false, /*.seqOffset=*/seqOffset};

                return std::make_tuple(std::move(bwt), std::move(csa));
            }();

            *this = KMerFMIndex{bwt, std::move(csa)};

        } else { // required 64bit SA required
            auto [bwt, csa] = [&]() {
                auto sa  = createSA64(inputText, threadNbr);
                auto bwt = createBWT64(inputText, sa);
                auto csa = TCSA{std::move(sa), samplingRate, inputSizes, /*.reverse=*/false, /*.seqOffset=*/seqOffset};

                return std::make_tuple(std::move(bwt), std::move(csa));
            }();

            *this = KMerFMIndex{bwt, std::move(csa)};
        }
    }
    auto operator=(KMerFMIndex const&) -> KMerFMIndex& = delete;
    auto operator=(KMerFMIndex&&) noexcept -> KMerFMIndex& = default;

    size_t size() const {
        return bwt.size();
    }

    auto locate(size_t idx) const -> std::tuple<ADEntry, size_t> {
        auto opt = csa.value(idx);
        size_t steps{};
        while(!opt) {
            auto symb = bwt.symbol(idx);
            idx = bwt.rank(idx, symb) + C[symb];
            steps += 1;
            opt = csa.value(idx);
        }
        return {*opt, steps};
    }

    auto single_locate_step(size_t idx) const -> std::optional<ADEntry> {
        return csa.value(idx);
    }

    template <typename Archive>
    void serialize(this auto&& self, Archive& ar) {
        ar(self.bwt, self.C, self.csa);
    }
};

}
