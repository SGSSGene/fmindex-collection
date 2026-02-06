// SPDX-FileCopyrightText: 2026 Simon Gene Gottlieb
// SPDX-License-Identifier: CC0-1.0
#include "allStrings.h"
#include "utils.h"

#include <fmindex-collection/string/PairedFlattenedBitvectors2L_b.h>

TEST_CASE("check if rank on strings with 'dual_limit' functions work", "[string][dual_limit]") {
    size_t const Sigma = 4;
    auto input1 = std::vector<uint8_t>{
        0, 1, 2, 3
    };
/*    auto input2 = std::vector<uint8_t>{};
    for (auto v1 : input1) {
        for (auto v2 : input1) {
            input2.push_back(Sigma*v1 + v2);
        }
    }*/

    auto input2 = generateText<0, 16>(1000);

    auto p = fmc::string::PairedFlattenedBitvectors2L_b<16, 2, 512, 65536>{input2};
//    auto p2 = fmc::string::PairedFlattenedBitvectors2L<16, 512, 65536>{input2};
//    (void)p2;
//    auto p = fmc::string::InterleavedBitvectorPrefix16<16>{input2};


    auto countR1 = [&](size_t idx, size_t symb) {
        size_t ct{};
        for (size_t i{0}; i < idx; ++i) {
            ct += (input2[i]/Sigma == symb);
        }
        return ct;
    };
    auto countPR1 = [&](size_t idx, size_t symb) {
        size_t ct{};
        for (size_t i{0}; i < idx; ++i) {
            ct += (input2[i]/Sigma < symb);
        }
        return ct;
    };

    auto countR2 = [&](size_t idx, size_t symb) {
        size_t ct{};
        for (size_t i{0}; i < idx; ++i) {
            ct += (input2[i] == symb);
        }
        return ct;
    };
    auto countPR2 = [&](size_t idx, size_t symb) {
        size_t ct{};
        for (size_t i{0}; i < idx; ++i) {
            ct += (input2[i] < symb);
        }
        return ct;
    };

/*    for (size_t i{0}; i <= input2.size(); ++i) {
        p.all_ranks_dual(i, input2.size(), [&](size_t symb, size_t rs1, size_t rs2, size_t prs1, size_t prs2) {
            CHECK(symb <= Sigma*Sigma);
            CHECK(rs1 == countR2(i, symb));
            CHECK(rs2 == countR2(input2.size(), symb));
            CHECK(prs1 == countPR2(i, symb));
            CHECK(prs2 == countPR2(input2.size(), symb));
        });
    }

    for (size_t i{0}; i <= input2.size(); ++i) {
        p.all_ranks_dual_limit<2>(i, input2.size(), [&](size_t symb, size_t rs1, size_t rs2, size_t prs1, size_t prs2) {
            CHECK(symb <= Sigma);
            CHECK(rs1 == countR1(i, symb));
            CHECK(rs2 == countR1(input2.size(), symb));
            CHECK(prs1 == countPR1(i, symb));
            CHECK(prs2 == countPR1(input2.size(), symb));
        });
    }*/

/*    for (size_t i{0}; i < input2.size(); ++i) {
        CHECK(p.symbol_limit<2>(i) == (input2[i] / Sigma));
    }*/
    for (size_t i{0}; i < input2.size(); ++i) {
        for (size_t s{0}; s < Sigma; ++s) {
            CHECK(p.rank_limit<2>(i, s) == countR1(i, s));
        }
    }
    for (size_t i{0}; i < input2.size(); ++i) {
        for (size_t s{0}; s <= Sigma; ++s) {
            CHECK(p.prefix_rank_limit<2>(i, s) == countPR1(i, s));
        }
    }

    for (size_t i{0}; i < input2.size(); ++i) {
        for (size_t s{0}; s < Sigma*Sigma; ++s) {
            auto [pr, r] = p.prefix_rank_and_rank(i, s);
            CHECK(pr == countPR2(i, s));
            CHECK(r == countR2(i, s));
        }
    }

    for (size_t i{0}; i < input2.size(); ++i) {
        for (size_t s{0}; s < Sigma; ++s) {
            auto [pr, r] = p.prefix_rank_and_rank_limit<2>(i, s);
            CHECK(pr == countPR1(i, s));
            CHECK(r == countR1(i, s));
        }
    }
}
