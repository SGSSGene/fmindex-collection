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
#include "CSA.h"

#include "utils/utils.h"

#include "search/SearchNg12.h"

#include "oss/generator/pigeon.h"
#include "oss/generator/h2.h"
#include "oss/generator/kucherov.h"
#include "oss/expand.h"

#include <fmt/format.h>

template <size_t Sigma, typename CSA, template <size_t> typename Table>
auto loadIndex(std::string path) {
    auto bwt       = readFile(path + ".bwt");
    auto bwtRev    = readFile(path + ".rev.bwt");
    auto csaBuffer = readFile(path + ".csa");

    CSA csa = [&]() {
        BitStack bitStack;
        std::vector<uint64_t> ssa;
        auto readLen = bitStack.read(csaBuffer.data(), csaBuffer.size());
        ssa.resize(bitStack.ones);
        memcpy(ssa.data(), csaBuffer.data() + readLen, ssa.size()* sizeof(uint64_t));
        return CSA(ssa, bitStack);
    }();
    auto index = BiFMIndex<Table<Sigma>>{bwt, bwtRev, csa};
    return index;
}

int main() {
    StopWatch watch;
    constexpr size_t Sigma = 5;

//    auto [bwt, bwtRev] = generateBWT<Sigma>(1ul<<20);
    auto index = loadIndex<Sigma, CSA, occtable::compact::OccTable>("/home/gene/hg38/short");

    std::vector<std::vector<uint8_t>> queries;
    {
        auto b = readFile("/home/gene/hg38/sampled_illumina_reads.txt");
//        auto b = readFile("/home/gene/hg38/short.queries");
        auto ptr = b.data();
        std::vector<uint8_t> query;
        while (ptr != (b.data() + b.size())) {
            if (*ptr == '\n') {
                queries.push_back(query);
                std::reverse(query.begin(), query.end());
                queries.push_back(query);
                query.clear();
            } else {
                if (*ptr == '$') query.push_back(0);
                else if (*ptr == 'A') query.push_back(1);
                else if (*ptr == 'C') query.push_back(2);
                else if (*ptr == 'G') query.push_back(3);
                else if (*ptr == 'T') query.push_back(4);
                else if (*ptr == 'N') query.push_back(5);
            }
            ++ptr;
        }
    }

    auto cursor = BiFMIndexCursor{index};
    std::cout << "start: " << cursor.lb << ", " << cursor.lbRev << " len: " << cursor.len << "\n";
    {
        StopWatch sw;
        size_t resultCt{};
        for (size_t i{0}; i < queries.size(); ++i) {
            auto const& q = queries[i];
            auto cursor = BiFMIndexCursor{index};
            for (size_t j{0}; j < q.size(); ++j) {
                cursor = cursor.extendRight(q[j]);
                if (cursor.empty()) {
                    break;
                }
            }
            if (!cursor.empty()) {
                resultCt += cursor.len;
            }
        }
        auto t = sw.reset();
        fmt::print("right: queries {}, took {}s, results: {}\n", queries.size(), t, resultCt);
    }

    {
        StopWatch sw;
        size_t resultCt{};
        for (size_t i{0}; i < queries.size(); ++i) {
            auto const& q = queries[i];
            auto cursor = BiFMIndexCursor{index};
            for (size_t j{0}; j < q.size(); ++j) {
                cursor = cursor.extendLeft(q[q.size() - j - 1]);
                if (cursor.empty()) {
                    break;
                }
            }
            if (!cursor.empty()) {
                resultCt += cursor.len;
            }
        }
        auto t = sw.reset();
        fmt::print("left: queries {}, took {}s, results: {}\n", queries.size(), t, resultCt);
    }



/*    fmt::print("one expansion\n");
    for (size_t i{1}; i < Sigma; ++i) {
        auto c2 = cursor.extendLeft(i);
        std::cout << i << " - start: " << c2.lb << "-" << c2.lb+c2.len << " or " << c2.lbRev << "-" << c2.lbRev+c2.len << " len: " << c2.len << "\n";
    }
    std::cout << "same as?\n";
    for (size_t i{1}; i < Sigma; ++i) {
        auto c2 = cursor.extendRight(i);
        std::cout << i << " - start: " << c2.lb << "-" << c2.lb+c2.len << " or " << c2.lbRev << "-" << c2.lbRev+c2.len << " len: " << c2.len << "\n";
    }

    fmt::print("two expansions\n");
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


    for (size_t k{0}; k<10; ++k)
    {
//        auto search_scheme = oss::expand(oss::generator::pigeon_trivial(0, k), queries[0].size());
        auto search_scheme = oss::expand(oss::generator::h2(k+2, 0, k), queries[0].size());
//        auto search_scheme = oss::expand(oss::generator::kucherov(k+2, 0, k), queries[0].size());
        for (size_t i{0}; i < search_scheme.size(); ++i) {
            auto& tree = search_scheme[i];
//            fmt::print("pi: ");
            for (size_t j{0}; j < tree.pi.size(); ++j) {
                tree.pi[j] -= 1;
//                fmt::print("{} ", tree.pi[j]);
            }
//            fmt::print("\n");
        }
        size_t resultCt{};
        StopWatch sw;
        std::vector<size_t> results{};
        search_ng12::search(index, queries, 0, search_scheme, [&](size_t idx, auto cursor) {
            for (size_t i{cursor.lb}; i < cursor.lb + cursor.len; ++i) {
                results.push_back(index.locate(i));
    //            fmt::print("found at {} (index: {})\n", index.locate(i), i);
//                fmt::print("found at {} (index: {})\n", csa.value(i).value_or(99999), i);
            }
            resultCt += cursor.len;
        });
        auto t = sw.reset();
        fmt::print("queries {}, took {}s, results: {}/{}\n", queries.size(), t, resultCt, results.size());
    }

    return 0;


/*
    auto time_bwtconstruction = watch.reset();
    std::cout << "bwt - construction time: "<< time_bwtconstruction << "s, length: " << index.size() << "\n";

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
    }*/

    return 0;
}


