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

#include "search/SearchPseudo.h"
#include "search/SearchNg12.h"
#include "search/SearchNg14.h"
#include "search/SearchNg15.h"
#include "search/SearchNg16.h"
#include "search/SearchNg17.h"
#include "search/SearchNg20.h"
#include "search/SearchNg21.h"

#include "oss/generator/pigeon.h"
#include "oss/generator/h2.h"
#include "oss/generator/kucherov.h"
#include "oss/generator/backtracking.h"
#include "oss/expand.h"

#include <fmt/format.h>
#include <unordered_set>

template <size_t Sigma, typename CB>
void visitAllTables(CB cb) {
    cb((occtable::naive::OccTable<Sigma>*)nullptr, "naive");
    cb((occtable::compact2::OccTable<Sigma>*)nullptr, "compact2");
    cb((occtable::compact2Aligned::OccTable<Sigma>*)nullptr, "compact2Aligned");
    cb((occtable::compactWavelet::OccTable<Sigma>*)nullptr, "compactWavelet");
    cb((occtable::sdsl_wt_bldc::OccTable<Sigma>*)nullptr, "sdsl_wt_bldc");
    cb((occtable::wavelet::OccTable<Sigma>*)nullptr, "wavelet");
    cb((occtable::compact::OccTable<Sigma>*)nullptr, "compact");
    cb((occtable::compactAligned::OccTable<Sigma>*)nullptr, "compactAligned");
    cb((occtable::compactPrefix::OccTable<Sigma>*)nullptr, "compactPrefix");
    cb((occtable::bitvector::OccTable<Sigma>*)nullptr, "bitvector");
    cb((occtable::bitvectorPrefix::OccTable<Sigma>*)nullptr, "bitvectorPrefix");
}


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

struct Query {
    std::string name;
    bool reverse;
};
template <size_t Sigma>
auto loadQueries(std::string path) {
    std::vector<std::vector<uint8_t>> queries;
    std::vector<Query> queryInfos;
    enum class Mode { Name, Sequence };
    {
        auto b = readFile(path);
        auto ptr = b.data();
        std::vector<uint8_t> query;
        auto mode = Mode::Name;
        if (*ptr != '>') {
            throw std::runtime_error("can't read fasta file");
        }
        while (ptr != (b.data() + b.size())) {
            if (mode == Mode::Name) {
                std::string name;
                if (ptr[0] != '>' or (ptr[1] != ' ')) {
                    throw std::runtime_error("expected '> '");
                }
                ++ptr;
                ++ptr;
                while (ptr != (b.data() + b.size()) and *ptr != '\n') {
                    name += *ptr;
                    ++ptr;
                }
                ++ptr;
                mode = Mode::Sequence;
                queryInfos.emplace_back(name, false);
                queryInfos.emplace_back(name, true);
            } else if (mode == Mode::Sequence) {
                if (*ptr == '>' || (ptr+1) == (b.data() + b.size())) {
                    queries.push_back(query);
                    std::reverse(query.begin(), query.end());
                    for (auto& c : query) {
                        if (c == 1) c = 4;
                        else if (c == 2) c = 3;
                        else if (c == 3) c = 2;
                        else if (c == 4) c = 1;
                    }
                    queries.push_back(query);
                    query.clear();
                    mode = Mode::Name;
                    if ((ptr+1) == (b.data() + b.size())) {
                        ++ptr;
                    }
                } else {
                    if (*ptr == '$') query.push_back(0);
                    else if (*ptr == 'A') query.push_back(1);
                    else if (*ptr == 'C') query.push_back(2);
                    else if (*ptr == 'G') query.push_back(3);
                    else if (*ptr == 'T') query.push_back(4);
                    else if (*ptr == 'N' and Sigma == 6) query.push_back(5);
                    else if (*ptr == '\n') {}
                    else {
                        throw std::runtime_error("unknown alphabet");
                    }
                    ++ptr;
                }
            }
        }
    }
    return std::make_tuple(queries, queryInfos);
}


int main() {
    StopWatch watch;
    constexpr size_t Sigma = 5;

     auto const [queries, queryInfos] = loadQueries<Sigma>("/home/gene/hg38/sampled_illumina_reads.fasta");
     fmt::print("loaded {} queries (incl reverse complements)\n", queries.size());

    fmt::print("{:15}: {:>10}  ({:>10} +{:>10} ) {:>10}    - results: {:>10}/{:>10}/{:>10}/{:>10} - mem: {:>13}\n", "name", "time_search + time_locate", "time_search", "time_locate", "(time_search+time_locate)/queries.size()", "resultCt", "results.size()", "uniqueResults.size()", "readIds.size()", "memory");


//    auto [bwt, bwtRev] = generateBWT<Sigma>(1ul<<20);
    visitAllTables<Sigma>([&]<template <size_t> typename Table>(Table<Sigma>*, std::string name) {
        auto mut_queries = queries;

        size_t s = Table<Sigma>::expectedMemoryUsage(3'000'000'000ul);
//        fmt::print("expected memory: {}\n", s);
        if (s > 1024*1024*1024*56ul) {
            fmt::print("{} skipping, to much memory\n", name);
            return;
        }
        auto index = loadIndex<Sigma, CSA, Table>("/home/gene/hg38/text.dna4");

        auto memory = index.memoryUsage();
        for (size_t k{0}; k<6; ++k)
        {
            if (k == 4 or k == 5) {
                mut_queries.resize(mut_queries.size() / 10);
            }

//            auto search_scheme = oss::expand(oss::generator::pigeon_opt(0, k), mut_queries[0].size());
//            auto search_scheme = oss::expand(oss::generator::pigeon_trivial(0, k), mut_queries[0].size());
            auto search_scheme = oss::expand(oss::generator::h2(k+2, 0, k), mut_queries[0].size());
//            auto search_scheme = oss::expand(oss::generator::kucherov(k+1, 0, k), mut_queries[0].size());
//            auto search_scheme = oss::expand(oss::generator::backtracking(1, 0, k), mut_queries[0].size());

    //
            for (size_t i{0}; i < search_scheme.size(); ++i) {
                auto& tree = search_scheme[i];
                for (size_t j{0}; j < tree.pi.size(); ++j) {
                    tree.pi[j] -= 1;
                }
            }
            size_t resultCt{};
            StopWatch sw;
            std::vector<std::tuple<size_t, size_t, size_t>> results{};
            std::vector<std::tuple<size_t, LeftBiFMIndexCursor<decltype(index)>, size_t>> resultCursors;

//            search_pseudo::search<true>(index, mut_queries, search_scheme, [&](size_t queryId, auto cursor, size_t errors) {
//            search_ng12::search(index, mut_queries, search_scheme, [&](size_t queryId, auto cursor, size_t errors) {
//            search_ng14::search(index, mut_queries, search_scheme, [&](size_t queryId, auto cursor, size_t errors) {
//            search_ng15::search(index, mut_queries, search_scheme, [&](size_t queryId, auto const& cursor, size_t errors) {
            search_ng16::search(index, mut_queries, search_scheme, [&](size_t queryId, auto const& cursor, size_t errors) {
//            search_ng20::search(index, mut_queries, search_scheme, [&](size_t queryId, auto cursor, size_t errors) {
//            search_ng21::search(index, mut_queries, search_scheme, [&](size_t queryId, auto cursor, size_t errors) {

                resultCursors.emplace_back(queryId, cursor, errors);
            });

            auto time_search = sw.reset();

            for (auto const& [queryId, cursor, e] : resultCursors) {
                for (size_t i{cursor.lb}; i < cursor.lb + cursor.len; ++i) {
                    results.emplace_back(queryId, index.locate(i), e);
                }
                resultCt += cursor.len;
            }
            auto time_locate = sw.reset();

            auto uniqueResults = [](auto list) {
                std::sort(begin(list), end(list));
                list.erase(std::unique(begin(list), end(list)), list.end());
                return list;
            }(results);
            std::unordered_set<size_t> readIds;
            for (auto [queryId, cursor, e] : resultCursors) {
                if (queryId > mut_queries.size()/2) {
                    readIds.insert(queryId - mut_queries.size() / 2);
                } else {
                    readIds.insert(queryId);
                }
            }

            fmt::print("{:15}: {:>10.3}s ({:>10.3}s+{:>10.3}s) {:>10.3}q/s - results: {:>10}/{:>10}/{:>10}/{:>10} - mem: {:>13}\n", name, time_search + time_locate, time_search, time_locate, (time_search+time_locate)/mut_queries.size(), resultCt, results.size(), uniqueResults.size(), readIds.size(), memory);
            {
                auto filename =fmt::format("out.k{}.ss{}.txt", k, name);
                /*auto ofs = fopen(filename.c_str(), "w");
                fmt::print(ofs, "identifier\tposition\tlength\tED\treverseComlpement\n");
                for (auto [queryId, pos, e] : results) {
                    auto const& qi = queryInfos[queryId];
                    fmt::print(ofs, "{}\t{}\t{}\t{}\t{}\n", qi.name, pos, 101, e, qi.reverse?1:0);
                }
                fclose(ofs);*/
            }
        }
    });
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


