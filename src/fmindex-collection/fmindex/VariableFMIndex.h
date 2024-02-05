// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../locate.h"
#include "../occtable/all.h"
#include "../rankvector/InterleavedBitvector.h"
#include "../search/SearchNoErrors.h"
#include "FMIndex.h"

#include <variant>

namespace fmindex_collection {

/**
 * Depending on the input it will choose a appropriate FMIndex
 */
struct VariableFMIndex {
    size_t Sigma{};

    std::array<uint8_t, 256> charToRankMapping{};

    template <size_t Sigma>
    using OccTable = occtable::Interleaved_16<Sigma>;

    using Index4  = FMIndex<OccTable<5>>;
    using Index5  = FMIndex<OccTable<6>>;
    using Index16 = FMIndex<OccTable<17>>;

    std::variant<std::monostate, Index4, Index5, Index16> index;

    VariableFMIndex() = default;

    VariableFMIndex(std::vector<std::string> const& _reference, size_t samplingRate, size_t threadNbr) {
        charToRankMapping.fill(255);

        // Scan and build up ranking
        // 1. set entries to 0 if they appear (branchless)
        for (auto const& str : _reference) {
            for (auto const& c : str) {
                charToRankMapping[static_cast<uint8_t>(c)] = 0;
            }
        }
        // 2. iterate over each entry and give them a uniq rank
        for (auto& c : charToRankMapping) {
            if (c == 0) {
                Sigma += 1;
                c = Sigma;
            }
        }

        // Convert _reference to compact rank representation
        auto reference = std::vector<std::vector<uint8_t>>{};
        reference.resize(_reference.size());
        for (size_t i{0}; i < _reference.size(); ++i) {
            auto const& str = _reference[i];
            auto& rankStr   = reference[i];

            rankStr.resize(str.size());
            for (size_t j{0}; j < str.size(); ++j) {
                rankStr[j] = charToRankMapping[str[j]];
            }
        }

        if (Sigma < 5) {
            index.emplace<Index4>(reference, samplingRate, threadNbr);
        } else if (Sigma < 6) {
            index.emplace<Index5>(reference, samplingRate, threadNbr);
        } else if (Sigma < 17) {
            index.emplace<Index16>(reference, samplingRate, threadNbr);
        } else {
            throw std::runtime_error{"No FMIndex available that can deal with more than 16 different characters"};
        }
    }

    auto search(std::string const& _query) const {
        // convert query to compact rank representation
        auto query = std::vector<uint8_t>{};
        query.resize(_query.size());
        for (size_t i{0}; i < _query.size(); ++i) {
            query[i] = charToRankMapping[_query[i]];
        }

        auto result = std::vector<std::tuple<size_t, size_t>>{};

        // check for invalid mapping characters
        {
            auto str = std::basic_string_view<uint8_t>{query.begin(), query.end()};
            auto pos = str.find(255);
            if (pos != std::string_view::npos) {
                // We will not find this query
                return result;
            }
        }

        std::visit([&]<typename I>(I const& index) {
            if constexpr (std::same_as<I, std::monostate>) {
                return;
            } else {
                auto cursor = search_no_errors::search(index, query);
                for (auto [seqId, pos] : LocateLinear{index, cursor}) {
                    result.emplace_back(seqId, pos);
                }
            }
        }, index);
        return result;
    }

    template <typename Archive>
    void save(Archive& ar) const {
        ar(size_t{1}); // Version 1
        ar(Sigma);
        ar(charToRankMapping);
        ar(size_t{index.index()});
        std::visit([&]<typename I>(I const& index) {
            if constexpr (std::same_as<I, std::monostate>) {
                return;
            } else {
                ar(index);
            }
        }, index);
    }

    template <typename Archive>
    void load(Archive& ar) {
        size_t version;
        ar(version);
        if (version == 1) {
            ar(Sigma);
            ar(charToRankMapping);
            size_t idx;
            ar(idx);
            if (idx == 1) index.emplace<Index4>();
            else if (idx == 2) index.emplace<Index5>();
            else if (idx == 3) index.emplace<Index16>();
            else throw std::runtime_error{"unknown index"};
            std::visit([&]<typename I>(I& index) {
                if constexpr (std::same_as<I, std::monostate>) {
                    return;
                } else {
                    ar(index);
                }
            }, index);
        } else {
            throw std::runtime_error{"unknown format"};
        }
    }
};

}
