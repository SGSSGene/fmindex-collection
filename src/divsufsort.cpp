#include "occtable/Bitvector.h"
#include "occtable/BitvectorPrefix.h"
#include "occtable/Compact.h"
#include "occtable/Compact2.h"
#include "occtable/Compact2Aligned.h"
#include "occtable/CompactAligned.h"
#include "occtable/CompactPrefix.h"
#include "occtable/CompactWavelet.h"
#include "occtable/Naive.h"
#include "occtable/Sdsl.h"
#include "occtable/Wavelet.h"

#include "concepts.h"
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
        std::cout << i << " " << int(bwt[i]) << "\t";
        for (size_t j{0}; j < table.Sigma; ++j) {
            std::cout << " " << table.rank(j, i);
        }
        std::cout << "\n";
    }
    std::cout << "\n\n";

    for (size_t i{0}; i < bwt.size(); ++i) {
        std::cout << i << " " << int(bwt[i]) << "\t";
        for (size_t j{0}; j < table.Sigma; ++j) {
            std::cout << " " << table.prefix_rank(j, i);
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
        xorshf96_reset();

        { // benchmark V1
            uint64_t a{};
            for (size_t i{0}; i < 10'000'000; ++i) {
                auto symb = xorshf96() % Table::Sigma;
                auto row = xorshf96() % table.size();
                a += table.rank(symb, row);
            }
            result.benchV1CheckSum = a;
            result.benchV1 = watch.reset();
        }
        { // benchmark V2
            size_t a = 0;
            for (size_t i{0}; i < 10'000'000; ++i) {
                auto symb = xorshf96() % Table::Sigma;
                auto row = xorshf96() % table.size();
                a += table.prefix_rank(symb, row);
            }
            result.benchV2CheckSum = a;
            result.benchV2 = watch.reset();
        }
        { // benchmark V3
            uint64_t jumps{1};
            uint64_t pos = table.rank(bwt[0], 0);
            uint64_t a{};
            while (pos != 0 && jumps/2 < bwt.size()) {
                jumps += 1;
                pos = table.rank(bwt[pos], pos);
                a += pos;
            }

            result.benchV3CheckSum = {jumps, a};
            result.benchV3 = watch.reset();
        }
        { // benchmark V4
            uint64_t a{};
            uint64_t jumps{1};
            uint64_t pos = table.rank(bwt[0], 0);
            while (pos != 0 && jumps/2 < bwt.size()) {
                jumps += 1;
                pos = table.rank(bwt[pos], pos);
                a += table.prefix_rank(bwt[pos], pos);
            }
            result.benchV4CheckSum = {jumps, a};
            result.benchV4 = watch.reset();
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
    text.resize(1ul<<30, '\0');
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
    }();*/
    auto bwt = readFile("/home/gene/hg38/text.dna5.bwt");
//    auto bwt = readFile("/home/gene/hg38/text.short.bwt");
//    bwt.resize(9);
/*    std::vector<uint8_t> bwt;
    bwt.resize(6, '\1');
    bwt[0] = 2;
    bwt[1] = 1;*/



    auto time_bwtconstruction = watch.reset();
    std::cout << "bwt - construction time: "<< time_bwtconstruction << "s, length: " << bwt.size() << "\n";

    auto results = std::vector<Result>{};
    using namespace occtable;
    results.emplace_back(benchmarkTable<naive::OccTable<Sigma>>("simple", bwt));
    results.emplace_back(benchmarkTable<compact::OccTable<Sigma>>("compact", bwt));
    results.emplace_back(benchmarkTable<compactAligned::OccTable<Sigma>>("compact_aligned", bwt));
    results.emplace_back(benchmarkTable<compactPrefix::OccTable<Sigma>>("compact_prefix", bwt));
    results.emplace_back(benchmarkTable<compact2::OccTable<Sigma>>("compact2", bwt));
    results.emplace_back(benchmarkTable<compact2Aligned::OccTable<Sigma>>("compact2_aligned", bwt));
    results.emplace_back(benchmarkTable<bitvector::OccTable<Sigma>>("bitvector", bwt));
    results.emplace_back(benchmarkTable<bitvectorPrefix::OccTable<Sigma>>("bitvectorPrefix", bwt));
    results.emplace_back(benchmarkTable<wavelet::OccTable<Sigma>>("bitvector_wavelet", bwt));
    results.emplace_back(benchmarkTable<compactWavelet::OccTable<Sigma>>("compact_wavelet2", bwt));
    results.emplace_back(benchmarkTable<occtable::sdsl::OccTable<Sigma>>("sdsl_wavelet", bwt));

    fmt::print(" {:^20} | {:^9} | {:^9} | {:^9} | {:^9} | {:^9} | {:^9} |"
               " {:^15} | {:^15} | {:^15} | {:^15} | {:^15} | {:^15} | {:^9}|\n",
                "name", "const", "mem(MB)", "bench 1", "bench 2", "bench 3", "bench 4",
                "check 1", "check 2", "check 3a", "check 3b", "check 4a", "check 4b", "total time");

    for (auto const& result : results) {
        fmt::print(" {:<20} | {:> 9.3} | {:> 9.0} | {:> 9.3} | {:> 9.3} | {:> 9.3} | {:> 9.3} |"
                   " {:>15} | {:>15} | {:>15} | {:>15} | {:>15} | {:>15} | {:> 9.3} |\n",
                   result.name,
                   result.constructionTime,
                   result.expectedMemory / 1024 / 1024,
                   result.benchV1,
                   result.benchV2,
                   result.benchV3,
                   result.benchV4,
                   result.benchV1CheckSum,
                   result.benchV2CheckSum,
                   result.benchV3CheckSum[0],
                   result.benchV3CheckSum[1],
                   result.benchV4CheckSum[0],
                   result.benchV4CheckSum[1],
                   result.totalTime);
    }

    return 0;
}


