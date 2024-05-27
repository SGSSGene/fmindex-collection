// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#include "../occtables/allTables.h"

#include <catch2/catch_all.hpp>
#include <fmindex-collection/fmindex/FMIndex.h>
#include <fmindex-collection/fmindex/merge.h>

TEST_CASE("checking merging of fmindices", "[FMIndex][merge]") {

    auto data1 = std::vector<std::vector<uint8_t>>{std::vector<uint8_t>{1, 1, 1, 1, 2, 2, 2, 2}};
    auto data2 = std::vector<std::vector<uint8_t>>{std::vector<uint8_t>{2, 1, 2, 1, 2, 1, 2, 2}};

    using OccTable = fmindex_collection::occtable::Bitvector<3>;
    using Index = fmindex_collection::FMIndex<OccTable>;

    auto index1 = Index{data1, /*.samplingRate =*/ 1, /*.threadNbr =*/ 1};
    auto index2 = Index{data2, /*.samplingRate =*/ 1, /*.threadNbr =*/ 1};

    std::cout << "idx: " << "\n";
    size_t idx1{};
    for (auto _ : data1[0]) {
        char c = index1.occ.symbol(idx1);
        idx1 = index1.occ.rank(idx1, c);
        std::cout << (int)c << " " << idx1 << "\n";
    }

    for (size_t i{0}; i < data1[0].size(); ++i) {
        auto idx0 = index1.occ.rank(i, 0);
        auto idx1 = index1.occ.rank(i, 1);
        auto idx2 = index1.occ.rank(i, 2);

        std::cout << idx0 << " " << idx1 << " " << idx2 << "\n";
    }

    auto index3 = merge(index1, index2);
    std::cout << "reconstruct:\n";
    auto texts = reconstructText(index3);
    for (auto const& t : texts) {
        for (auto c : t) {
            std::cout << (int)c;
        }
        std::cout << "\n";
    }
}
