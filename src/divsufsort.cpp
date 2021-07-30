#include "occtable/Bitvector.h"
#include "occtable/BitvectorPrefix.h"
#include "occtable/Compact.h"
#include "occtable/Compact2.h"
#include "occtable/Compact2Aligned.h"
#include "occtable/CompactAligned.h"
#include "occtable/CompactPrefix.h"
#include "occtable/CompactWavelet.h"
#include "occtable/Naive.h"
#include "occtable/Sdsl_wt_bldc.h"
#include "occtable/Wavelet.h"

#include "FMIndex.h"
#include "BiFMIndex.h"

#include "random.h"
#include "StopWatch.h"

#include <divsufsort64.h>

#include <array>
#include <bitset>
#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <fmt/format.h>


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


template <OccTable Table, typename T>
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
template <OccTable Table, typename T>
auto benchmarkTable(std::string name, T const& bwt) -> Result {
    StopWatch allTime;

    Result result;
    result.name = name;
    size_t s = Table::expectedMemoryUsage(bwt.size());

    result.expectedMemory = s;

    if (s < 1024*1024*1024*8ul) {
        StopWatch watch;
        auto table = Table{bwt};

        result.constructionTime = watch.reset();
        printOccTable(table, bwt);
        { // benchmark V1
            xorshf96_reset();
            uint64_t a{};
            for (size_t i{0}; i < 10'000'000; ++i) {
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
            for (size_t i{0}; i < 10'000'000; ++i) {
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

int main() {
    StopWatch watch;
    constexpr size_t Sigma = 6;

/*    std::string text;
    text.resize(1ul<<28, '\0');
    for (size_t i{0}; i < text.size(); ++i) {
        text[i] = (xorshf96() % (Sigma-1))+1;
    }
    text.back() = '\0';

    auto time_generation = watch.reset();
    std::cout << "text size: " << text.size()/1024/1024 << " million chars, "<<  text.size()/1024/1024*std::log(Sigma)/std::log(2) << " million bits\n";
    std::cout << "text generation: "<< time_generation << "s\n";

    auto bwt = [&]() {
        std::vector<int64_t> sa;
        sa.resize(text.size());
        auto error = divsufsort64((uint8_t const*)text.data(), sa.data(), text.size());
        if (error != 0) {
            throw std::runtime_error("some error while creating the suffix array");
        }

        auto time_saconstruction = watch.reset();
        std::cout << "sa - construction time: "<< time_saconstruction << "s\n";

        return construct_bwt_from_sa(sa, text);
    }();
    auto bwtRev = [&]() {
        std::vector<int64_t> sa;
        std::reverse(text.begin(), text.end());
        sa.resize(text.size());
        auto error = divsufsort64((uint8_t const*)text.data(), sa.data(), text.size());
        if (error != 0) {
            throw std::runtime_error("some error while creating the suffix array");
        }

        auto time_saconstruction = watch.reset();
        std::cout << "sa - rev construction time: "<< time_saconstruction << "s\n";

        return construct_bwt_from_sa(sa, text);
    }();*/

    auto bwt = readFile("/home/gene/hg38/text.dna5.bwt");
    auto bwtRev = readFile("/home/gene/hg38/text.dna5.rev.bwt");

//    auto bwt = readFile("/home/gene/hg38/text.short.bwt");
//    bwt.resize(9);
/*    std::vector<uint8_t> bwt;
    bwt.resize(6, '\1');
    bwt[0] = 2;
    bwt[1] = 1;
    bwt[2] = 3;
    bwt[3] = 3;*/

/*    auto index = BiFMIndex<occtable::compact::OccTable<Sigma>>{bwt, bwtRev};
    auto cursor = BiFMIndexCursor{index};
    std::cout << "start: " << cursor.lb << ", " << cursor.lbRev << " len: " << cursor.len << "\n";
    for (size_t i{1}; i < Sigma; ++i) {
        for (size_t j{1}; j < Sigma; ++j) {
            auto c = cursor.extendLeft(j);
            auto c2 = c.extendLeft(i);
            std::cout << j<<" " <<i << " - start: " << c2.lb << "-" << c2.lb+c2.len << " or " << c2.lbRev << "-" << c2.lbRev+c2.len << " len: " << c2.len << "\n";
        }
    }
    std::cout << "same as?\n";
    for (size_t i{1}; i < Sigma; ++i) {
        for (size_t j{1}; j < Sigma; ++j) {
            auto c = cursor.extendRight(j);
            auto c2 = c.extendLeft(i);
            std::cout << j<<" " <<i << " - start: " << c2.lb << "-" << c2.lb+c2.len << " or " << c2.lbRev << "-" << c2.lbRev+c2.len << " len: " << c2.len << "\n";
        }
    }
    std::cout << "same as?\n";
    for (size_t i{1}; i < Sigma; ++i) {
        for (size_t j{1}; j < Sigma; ++j) {
            auto c = cursor.extendLeft(i);
            auto c2 = c.extendRight(j);
            std::cout << j<<" " <<i << " - start: " << c2.lb << "-" << c2.lb+c2.len << " or " << c2.lbRev << "-" << c2.lbRev+c2.len << " len: " << c2.len << "\n";
        }
    }
    std::cout << "same as?\n";
    for (size_t i{1}; i < Sigma; ++i) {
        for (size_t j{1}; j < Sigma; ++j) {
            auto c = cursor.extendRight(i);
            auto c2 = c.extendRight(j);
            std::cout << j<<" " <<i << " - start: " << c2.lb << "-" << c2.lb+c2.len << " or " << c2.lbRev << "-" << c2.lbRev+c2.len << " len: " << c2.len << "\n";
        }
    }*/


//    return 0;



    auto time_bwtconstruction = watch.reset();
    std::cout << "bwt - construction time: "<< time_bwtconstruction << "s, length: " << bwt.size() << "\n";

    auto results = std::vector<Result>{};
    using namespace occtable;
    results.emplace_back(benchmarkTable<naive::OccTable<Sigma>>("naive", bwt));
    results.emplace_back(benchmarkTable<compact::OccTable<Sigma>>("compact", bwt));
    results.emplace_back(benchmarkTable<compactAligned::OccTable<Sigma>>("compactAligned", bwt));
    results.emplace_back(benchmarkTable<compactPrefix::OccTable<Sigma>>("compactPrefix", bwt));
    results.emplace_back(benchmarkTable<compact2::OccTable<Sigma>>("compact2", bwt));
    results.emplace_back(benchmarkTable<compact2Aligned::OccTable<Sigma>>("compact2Aligned", bwt));
    results.emplace_back(benchmarkTable<bitvector::OccTable<Sigma>>("bitvector", bwt));
    results.emplace_back(benchmarkTable<bitvectorPrefix::OccTable<Sigma>>("bitvectorPrefix", bwt));
    results.emplace_back(benchmarkTable<wavelet::OccTable<Sigma>>("wavelet", bwt));
    results.emplace_back(benchmarkTable<compactWavelet::OccTable<Sigma>>("compactWavelet", bwt));
    results.emplace_back(benchmarkTable<sdsl_wt_bldc::OccTable<Sigma>>("sdsl_wt_bldc", bwt));

    fmt::print(" {:^20} | {:^9} | {:^9} | {:^9} | {:^9} | {:^9} | {:^9} | {:^9} | {:^9} |"
               " {:^15} | {:^15} | {:^15} | {:^15} | {:^15} | {:^15} | {:^15} | {:^15} | {:^9} | \n",
                "name", "const", "mem(MB)", "bench 1", "bench 2", "bench 3", "bench 4", "bench 5", "bench 6",
                "check 1", "check 2", "check 3a", "check 3b", "check 4a", "check 4b", "check 5", "check 6a", "check 6b", "total time");

    for (auto const& result : results) {
        fmt::print(" {:<20} | {:> 9.3} | {:> 9.0} | {:> 9.3} | {:> 9.3} | {:> 9.3} | {:> 9.3} | {:> 9.3} | {:> 9.3} |"
                   " {:>15} | {:>15} | {:>15} | {:>15} | {:>15} | {:>15} | {:>15} | {:>15} | {:>15} | {:> 9.3} |\n",
                   result.name,
                   result.constructionTime,
                   result.expectedMemory / 1024 / 1024,
                   result.benchV1,
                   result.benchV2,
                   result.benchV3,
                   result.benchV4,
                   result.benchV5,
                   result.benchV6,
                   result.benchV1CheckSum,
                   result.benchV2CheckSum,
                   result.benchV3CheckSum[0],
                   result.benchV3CheckSum[1],
                   result.benchV4CheckSum[0],
                   result.benchV4CheckSum[1],
                   result.benchV5CheckSum,
                   result.benchV6CheckSum[0],
                   result.benchV6CheckSum[1],
                   result.totalTime);
    }

    return 0;
}


