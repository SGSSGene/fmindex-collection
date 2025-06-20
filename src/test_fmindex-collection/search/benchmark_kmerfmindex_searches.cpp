// SPDX-FileCopyrightText: 2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <catch2/catch_all.hpp>
#include <fmindex-collection/fmindex/KMerFMIndex.h>
#include <fmindex-collection/string/all.h>
#include <fmindex-collection/search/KMerSearch.h>
#include <nanobench.h>

TEST_CASE("benchmark searches with errors", "[searches][!benchmark][kmerfmindex]") {
    SECTION("benchmarking") {
        using String = fmc::string::InterleavedBitvector16<256>;
        using Index = fmc::KMerFMIndex<String, 8>;

        srand(0);

        auto generateSequence = [](size_t l) {
            auto seq = std::vector<uint8_t>{};
            seq.reserve(l);
            for (size_t i{0}; i < l; ++i) {
                seq.emplace_back(1 + rand()%4);
            }
            return seq;
        };

        auto generateSequences = [generateSequence](size_t n, size_t l) {
            auto ref = std::vector<std::vector<uint8_t>>{};
            ref.reserve(n);
            for (size_t i{0}; i < n; ++i) {
                ref.emplace_back(generateSequence(l));
            }
            return ref;
        };

        auto generateRead = [](std::vector<std::vector<uint8_t>> const& ref, size_t l, size_t e) {
            auto subjId = rand() % ref.size();
            auto& r = ref[subjId];

            auto edit = std::vector<char>{};
            edit.resize(l, 'M');
            while (e > 0) {
                auto t = rand() % 3;
                if (t == 0) { // substitution
                    auto pos2 = rand() % edit.size();
                    if (edit[pos2] != 'M') continue;
                    edit[pos2] = 'S';
                    --e;
                } else if (t == 1) { // insertion
                    auto pos2 = rand() % edit.size();
                    if (edit[pos2] != 'M') continue;
                    edit[pos2] = 'I';
                    --e;
                } else if (t == 2) { // deletion
                    auto pos2 = rand() % (edit.size()+1);
                    edit.insert(edit.begin() + pos2, 'D');
                    --e;
                }
            }

            auto pos = rand() % (r.size()-l-e);
            auto read = std::vector<uint8_t>{};
            for (auto t : edit) {
                if (t == 'M') {
                    read.push_back(r[pos]);
                    ++pos;
                } else if (t == 'S') {
                    auto c = ((r[pos]-1) + (rand() % 3)) % 4 + 1;
                    read.push_back(c);
                    ++pos;
                } else if (t == 'I') {
                    auto c = (rand() % 4) + 1;
                    read.push_back(c);
                } else if (t == 'D') {
                    ++pos;
                }
            }
            auto gt = std::make_tuple(subjId, pos);
            return std::make_tuple(gt, read);
        };

        auto generateReads = [generateRead](size_t n, std::vector<std::vector<uint8_t>> const& ref, size_t l, size_t e) {
            auto reads = std::vector<std::vector<uint8_t>>{};
            auto expected = std::vector<std::tuple<size_t, size_t, size_t>>{};
            reads.reserve(n);
            for (size_t i{0}; i < n; ++i) {
                auto [gt, read] = generateRead(ref, l, e);
                reads.emplace_back(read);
                expected.emplace_back(i, std::get<0>(gt), std::get<1>(gt));
            }
            return std::make_tuple(expected, reads);
        };

        // generate reference
        auto ref = generateSequences(100, 1'000'000);

        // generate reads
        for (size_t errors{0}; errors < 1; ++errors) {
            size_t len = 150;
            auto [expected, reads] = generateReads(1000, ref, len, errors);

            auto index = Index{ref, /*samplingRate*/1, /*threadNbr*/1};
            static auto bench = ankerl::nanobench::Bench();
            bench.batch(reads.size())
                 .relative(true);

            {
                bench.run("search kmersearch - error " + std::to_string(errors), [&]() {
                    fmc::kmersearch::search(index, reads, [&](auto qidx, auto cursor) {
                        (void)errors;
                        (void)qidx;
                        ankerl::nanobench::doNotOptimizeAway(cursor);
                    });
                });
            }
        }
    }
}
