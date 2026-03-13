// SPDX-FileCopyrightText: 2026 Simon Gene Gottlieb
// SPDX-License-Identifier: CC0-1.0
#include "allStrings.h"
#include "utils.h"

#include <fmindex-collection/string/PairedFlattenedBitvectors2LPartialSymb.h>

TEST_CASE("check if rank on strings with 'dual_limit' functions work", "[string][dual_limit]") {
    size_t const Sigma = 4;
    auto input1 = std::vector<uint8_t>{
        0, 1, 2, 3
    };

    auto input2 = generateText<0, 64>(1000);

    auto p = fmc::string::PairedFlattenedBitvectors2LPartialSymb<64, 2, 512, 65536>{input2};
//
    auto countR1_arr  = std::vector<std::array<size_t, 4>>{};
    {
        auto ct = std::array<size_t, 4>{};
        for (size_t i{0}; i < input2.size(); ++i) {
            countR1_arr.push_back(ct);
            ct[input2[i]/Sigma/Sigma] += 1;
        }
        countR1_arr.push_back(ct);
    }
    auto countR2_arr  = std::vector<std::array<size_t, 16>>{};
    {
        auto ct = std::array<size_t, 16>{};
        for (size_t i{0}; i < input2.size(); ++i) {
            countR2_arr.push_back(ct);
            ct[input2[i]/Sigma] += 1;
        }
        countR2_arr.push_back(ct);
    }

    auto countR3_arr  = std::vector<std::array<size_t, 64>>{};
    {
        auto ct = std::array<size_t, 64>{};
        for (size_t i{0}; i < input2.size(); ++i) {
            countR3_arr.push_back(ct);
            ct[input2[i]] += 1;
        }
        countR3_arr.push_back(ct);
    }

    auto countR1 = [&](size_t idx, size_t symb) {
        return countR1_arr[idx][symb];
    };
    auto countPR1 = [&](size_t idx, size_t symb) {
        size_t ct{};
        for (size_t i{0}; i < symb; ++i) {
            ct += countR1(idx, i);
        }
        return ct;
    };

    auto countR2 = [&](size_t idx, size_t symb) {
        return countR2_arr[idx][symb];
    };
    auto countPR2 = [&](size_t idx, size_t symb) {
        size_t ct{};
        for (size_t i{0}; i < symb; ++i) {
            ct += countR2(idx, i);
        }
        return ct;
    };

    auto countR3 = [&](size_t idx, size_t symb) {
        return countR3_arr[idx][symb];
    };
    auto countPR3 = [&](size_t idx, size_t symb) {
        size_t ct{};
        for (size_t i{0}; i < symb; ++i) {
            ct += countR3(idx, i);
        }
        return ct;
    };

    for (size_t i{0}; i <= input2.size(); ++i) {
        p.all_ranks_dual_limit<6>(i, input2.size(), [&](size_t symb, size_t rs1, size_t rs2, size_t prs1, size_t prs2) {
            CHECK(symb < Sigma*Sigma*Sigma);
            CHECK(rs1 == countR3(i, symb));
            CHECK(rs2 == countR3(input2.size(), symb));
            CHECK(prs1 == countPR3(i, symb));
            CHECK(prs2 == countPR3(input2.size(), symb));
        });
    }

    for (size_t i{0}; i <= input2.size(); ++i) {
        p.all_ranks_dual_limit<4>(i, input2.size(), [&](size_t symb, size_t rs1, size_t rs2, size_t prs1, size_t prs2) {
            CHECK(symb < Sigma*Sigma);
            CHECK(rs1 == countR2(i, symb));
            CHECK(rs2 == countR2(input2.size(), symb));
            CHECK(prs1 == countPR2(i, symb));
            CHECK(prs2 == countPR2(input2.size(), symb));
        });
    }

    for (size_t i{0}; i <= input2.size(); ++i) {
        p.all_ranks_dual_limit<2>(i, input2.size(), [&](size_t symb, size_t rs1, size_t rs2, size_t prs1, size_t prs2) {
            CHECK(symb < Sigma);
            CHECK(rs1 == countR1(i, symb));
            CHECK(rs2 == countR1(input2.size(), symb));
            CHECK(prs1 == countPR1(i, symb));
            CHECK(prs2 == countPR1(input2.size(), symb));
        });
    }

    for (size_t i{0}; i < input2.size(); ++i) {
        CHECK(p.symbol(i) == input2[i]);
        CHECK(p.symbol_limit<2>(i) == input2[i]/Sigma/Sigma);
        CHECK(p.symbol_limit<4>(i) == input2[i]/Sigma);
        CHECK(p.symbol_limit<6>(i) == input2[i]);
    }

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
            CHECK(p.rank_limit<4>(i, s) == countR2(i, s));
        }
    }
    for (size_t i{0}; i < input2.size(); ++i) {
        for (size_t s{0}; s <= Sigma*Sigma; ++s) {
            CHECK(p.prefix_rank_limit<4>(i, s) == countPR2(i, s));
        }
    }


    for (size_t i{0}; i < input2.size(); ++i) {
        for (size_t s{0}; s < Sigma*Sigma*Sigma; ++s) {
            auto [pr, r] = p.prefix_rank_and_rank(i, s);
            CHECK(pr == countPR3(i, s));
            CHECK(r == countR3(i, s));
        }
    }

    for (size_t i{0}; i < input2.size(); ++i) {
        for (size_t s{0}; s < Sigma; ++s) {
            {
                auto [pr, r] = p.prefix_rank_and_rank_limit<2>(i, s);
                CHECK(pr == countPR1(i, s));
                CHECK(r == countR1(i, s));
            }
            {
                auto [pr, r] = p.prefix_rank_and_rank_limit<4>(i, s);
                CHECK(pr == countPR2(i, s));
                CHECK(r == countR2(i, s));
            }
        }
    }
}
