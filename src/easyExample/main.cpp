// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#include <fmindex-collection/fmindex-collection.h>
#include <fmindex-collection/search/all.h>
#include <fmindex-collection/search_scheme/expand.h>
#include <fmindex-collection/search_scheme/generator/all.h>

#include <cereal/archives/binary.hpp>
#include <filesystem>
#include <fmt/format.h>
#include <fstream>
#include <unordered_set>

using namespace fmc;

constexpr size_t Sigma = 5;

template <size_t Sigma>
using String = string::InterleavedBitvector16<Sigma>;

template <typename Index>
void saveIndex(Index const& _index, std::filesystem::path _fileName) {
    auto ofs     = std::ofstream(_fileName, std::ios::binary);
    auto archive = cereal::BinaryOutputArchive{ofs};
    archive(_index);
}

template <typename Index>
auto loadIndex(std::filesystem::path _fileName) {
    auto ifs     = std::ifstream(_fileName, std::ios::binary);
    auto archive = cereal::BinaryInputArchive{ifs};
    auto index = Index{};
    archive(index);
    return index;
}


int main(int argc, char const* const* argv) {
    (void)argc;
    (void)argv;
    auto reference = std::vector<std::vector<uint8_t>>{{1, 3, 1, 4}, {2, 1, 4, 2, 3}};
    {
        std::cout << "\nBiFMIndex:\n";
        auto index = FMIndex<Sigma>{reference, /*samplingRate*/16, /*threadNbr*/1};
        auto queries = std::vector<std::vector<uint8_t>>{{1, 3}, {4, 2}};

        search_backtracking::search(index, queries, 0, [&](size_t queryId, auto cursor, size_t errors) {
            (void)errors;
            fmt::print("found something {} {}\n", queryId, cursor.count());
            for (auto i : cursor) {
                auto [entry, offset] = index.locate(i);
                auto [chr, pos] = entry;
                fmt::print("chr/pos: {}/{}\n", chr, pos+offset);
            }
        });
    }
    {
        std::cout << "\nFMIndex:\n";
        auto index = FMIndex<Sigma>{reference, /*samplingRate*/16, /*threadNbr*/1};
        auto queries = std::vector<std::vector<uint8_t>>{{1, 3}, {4, 2}};

        search_backtracking::search(index, queries, 0, [&](size_t queryId, auto cursor, size_t errors) {
            (void)errors;
            fmt::print("found something {} {}\n", queryId, cursor.count());
            for (auto i : cursor) {
                auto [entry, offset] = index.locate(i);
                auto [chr, pos] = entry;
                fmt::print("chr/pos: {}/{}\n", chr, pos+offset);
            }
        });
    }

    {
        std::cout << "\nReverseFMIndex:\n";
        auto index = ReverseFMIndex<String<Sigma>>{reference, /*samplingRate*/16, /*threadNbr*/1};
        auto queries = std::vector<std::vector<uint8_t>>{{3, 1}, {2, 4}};
        search_backtracking::search(index, queries, 0, [&](size_t queryId, auto cursor, size_t errors) {
            (void)errors;
            fmt::print("found something {} {}\n", queryId, cursor.count());
            for (auto i : cursor) {
                auto [entry, offset] = index.locate(i);
                auto [chr, pos] = entry;
                pos = pos+offset -  queries[queryId].size();
                fmt::print("chr/pos: {}/{}\n", chr, pos);
            }
        });

    }



/*    auto search_scheme = search_scheme::expand(search_scheme::generator::pigeon_opt(0, 0), queries[0].size());

    search_pseudo::search<false>(index, queries, search_scheme, [&](size_t queryId, auto cursor, size_t errors) {
        (void)errors;
        fmt::print("found something {} {}\n", queryId, cursor.count());
        for (auto i{begin(cursor)}; i < end(cursor); ++i) {
            auto [chr, pos] = index.locate(i);
            fmt::print("chr/pos: {}/{}\n", chr, pos);
        }
    });*/

    return 0;
}


