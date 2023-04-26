#pragma once

#include "StopWatch.h"
#include "random.h"

#include <fmindex-collection/occtable/concepts.h>
#include <fmindex-collection/utils.h>

#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>



inline auto construct_bwt_from_sa(std::vector<int64_t> const& sa, std::string_view const& text) -> std::vector<uint8_t> {
    assert(sa.size() == text.size());
    std::vector<uint8_t> bwt;
    bwt.resize(text.size());
    for (size_t i{0}; i < sa.size(); ++i) {
        bwt[i] = text[(sa[i] + text.size()- 1) % text.size()];
    }
    return bwt;
}
inline auto readFile(std::filesystem::path const& file) -> std::vector<uint8_t> {
    auto ifs = std::ifstream{file, std::ios::binary};
    ifs.seekg(0, std::ios::end);
    std::size_t size = ifs.tellg();
    auto buffer = std::vector<uint8_t>(size);
    ifs.seekg(0, std::ios::beg);
    ifs.read(reinterpret_cast<char*>(buffer.data()), buffer.size());
    return buffer;
}

template <size_t Sigma>
auto generateBWT(size_t len) -> std::tuple<std::vector<uint8_t>, std::vector<uint8_t>> {
    StopWatch watch;
    std::string text;
    text.resize(len, '\0');
    for (size_t i{0}; i < text.size(); ++i) {
        text[i] = (xorshf96() % (Sigma-1))+1;
    }
    text.back() = '\0';

    auto time_generation = watch.reset();
    std::cout << "text size: " << text.size()/1024/1024 << " million chars, "<<  text.size()/1024/1024*std::log(Sigma)/std::log(2) << " million bits\n";
    std::cout << "text generation: "<< time_generation << "s\n";

    auto bwt = [&]() {
        auto sa = fmindex_collection::createSA({reinterpret_cast<uint8_t const*>(text.data()), text.size()}, 1);
        auto time_saconstruction = watch.reset();
        std::cout << "sa - construction time: "<< time_saconstruction << "s\n";

        return construct_bwt_from_sa(sa, text);
    }();
    auto bwtRev = [&]() {
        std::ranges::reverse(text);
        auto sa = fmindex_collection::createSA({reinterpret_cast<uint8_t const*>(text.data()), text.size()}, 1);
        auto time_saconstruction = watch.reset();
        std::cout << "sa - rev construction time: "<< time_saconstruction << "s\n";

        return construct_bwt_from_sa(sa, text);
    }();
    return {bwt, bwtRev};
}


template <fmindex_collection::OccTable Table, typename T>
void printOccTable(Table const& table, T const& bwt) {
    if (bwt.size() > 128) return;
    for (size_t i{0}; i < bwt.size(); ++i) {
        auto r = std::get<0>(table.all_ranks(i));
        std::cout << i << " " << int(bwt[i]) << "\t";
        for (size_t j{0}; j < table.Sigma; ++j) {
            std::cout << " " << r[j];
//            std::cout << " " << table.rank(i, j);
        }
        std::cout << "\n";
    }
    std::cout << "\n\n";

    for (size_t i{0}; i < bwt.size(); ++i) {
        std::cout << i << " " << int(bwt[i]) << "\t";
        for (size_t j{0}; j < table.Sigma; ++j) {
            std::cout << " " << std::get<1>(table.all_ranks(i))[j];
        }
        std::cout << "\n";
    }
}

struct Result {
    std::string name;
    double expectedMemory                  = std::numeric_limits<double>::quiet_NaN();
    double constructionTime                = std::numeric_limits<double>::quiet_NaN();
    double benchV1                         = std::numeric_limits<double>::quiet_NaN();
    size_t benchV1CheckSum                 = 0;
    double benchV2                         = std::numeric_limits<double>::quiet_NaN();
    size_t benchV2CheckSum                 = 0;
    double benchV3                         = std::numeric_limits<double>::quiet_NaN();
    std::array<size_t, 2> benchV3CheckSum  = {0, 0};
    double benchV4                         = std::numeric_limits<double>::quiet_NaN();
    std::array<size_t, 2> benchV4CheckSum  = {0, 0};
    double benchV5                         = std::numeric_limits<double>::quiet_NaN();
    size_t benchV5CheckSum                 = 0;
    double benchV6                         = std::numeric_limits<double>::quiet_NaN();
    std::array<size_t, 2> benchV6CheckSum  = {0, 0};
    double totalTime                       = std::numeric_limits<double>::quiet_NaN();
};
template <fmindex_collection::OccTable Table, typename T>
auto benchmarkTable(std::string name, T const& bwt) -> Result {
    StopWatch allTime;

    Result result;
    result.name = name;
    size_t s = Table::expectedMemoryUsage(bwt.size());

    result.expectedMemory = s;

    if (s < 1024ull*1024*1024*8ull) {
        StopWatch watch;
        auto table = Table{bwt};

        result.constructionTime = watch.reset();
        printOccTable(table, bwt);
        { // benchmark V1
            xorshf96_reset();
            uint64_t a{};
            for (size_t i{0}; i < 10'000'000ull; ++i) {
                auto symb = xorshf96() % Table::Sigma;
                auto row = xorshf96() % table.size();
                a += table.rank(row, symb);
            }
            result.benchV1CheckSum = a;
            result.benchV1 = watch.reset();
        }
        { // benchmark V2
            xorshf96_reset();
            size_t a = 0;
            for (size_t i{0}; i < 10'000'000ull; ++i) {
                auto symb = xorshf96() % Table::Sigma;
                auto row = xorshf96() % table.size();
                a += table.prefix_rank(row, symb);
            }
            result.benchV2CheckSum = a;
            result.benchV2 = watch.reset();
        }
        { // benchmark V3
            uint64_t jumps{1};
            uint64_t pos = table.rank(0, bwt[0]);
            uint64_t a{};
            while (pos != 0 && jumps/2 < bwt.size()) {
                jumps += 1;
                pos = table.rank(pos, bwt[pos]);
                a += pos;
            }

            result.benchV3CheckSum = {jumps, a};
            result.benchV3 = watch.reset();
        }
        { // benchmark V4
            uint64_t a{};
            uint64_t jumps{1};
            uint64_t pos = table.rank(0, bwt[0]);
            while (pos != 0 && jumps/2 < bwt.size()) {
                jumps += 1;
                a += table.prefix_rank(pos, bwt[pos]);
                pos = table.rank(pos, bwt[pos]);
            }
            result.benchV4CheckSum = {jumps, a};
            result.benchV4 = watch.reset();
        }
        { // benchmark V5
            xorshf96_reset();
            uint64_t a{};
            for (size_t i{0}; i < 10'000'000; ++i) {
                auto symb = xorshf96() % Table::Sigma;
                auto row = xorshf96() % table.size();
                auto [rs, prs] = table.all_ranks(row);
                for (size_t j{0}; j < Table::Sigma; ++j) {
                    a += rs[j] + prs[j];
                }
            }
            result.benchV5CheckSum = a;
            result.benchV5 = watch.reset();
        }
        { // benchmark V6
            uint64_t jumps{1};
            uint64_t pos = table.rank(0, bwt[0]);
            uint64_t a{};
            while (pos != 0 && jumps/2 < bwt.size()) {
                jumps += 1;
                auto [rs, prs] = table.all_ranks(pos);
                pos = rs[bwt[pos]];

                for (size_t j{0}; j < Table::Sigma; ++j) {
                    a += rs[j] + prs[j];
                }
            }

            result.benchV6CheckSum = {jumps, a};
            result.benchV6 = watch.reset();
        }



    }
    result.totalTime = allTime.reset();
    std::cout << "finished " << result.name << "\n";
    return result;
}

