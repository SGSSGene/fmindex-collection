#include "utils/utils.h"

#include <fmindex-collection/fmindex-collection.h>
#include <fmindex-collection/occtable/all.h>
#include <fmindex-collection/search/all.h>
#include <fmindex-collection/locate.h>
#include <search_schemes/generator/all.h>
#include <search_schemes/expand.h>

#include <cereal/archives/binary.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/vector.hpp>

#include <fmt/format.h>
#include <unordered_set>

using namespace fmindex_collection;

template <size_t Sigma, typename CB>
void visitAllTables(CB cb) {
//    cb(std::type_identity<occtable::naive::OccTable<Sigma>>{}, "naive");
//    cb(std::type_identity<occtable::bitvector::OccTable<Sigma>>{}, "bitvector");
//    cb(std::type_identity<occtable::interleaved16::OccTable<Sigma>>{}, "interleaved16");
//    cb(std::type_identity<occtable::interleaved32::OccTable<Sigma>>{}, "interleaved32");
    cb(std::type_identity<occtable::interleavedEPR16::OccTable<Sigma>>{}, "interleavedEPR16");
    cb(std::type_identity<occtable::interleavedEPR32::OccTable<Sigma>>{}, "interleavedEPR32");
    cb(std::type_identity<occtable::interleavedEPR16V2::OccTable<Sigma>>{}, "interleavedEPR16V2");
    cb(std::type_identity<occtable::interleavedEPR32V2::OccTable<Sigma>>{}, "interleavedEPR32V2");
    cb(std::type_identity<occtable::eprV4::OccTable<Sigma>>{}, "EPRV4");
    cb(std::type_identity<occtable::eprV5::OccTable<Sigma>>{}, "EPRV5");
//    cb(std::type_identity<occtable::interleavedWavelet::OccTable<Sigma>>{}, "compactWavelet");
//    cb(std::type_identity<occtable::compactWaveletAligned::OccTable<Sigma>>{}, "compactWaveletAligned");
//    cb(std::type_identity<occtable::compact2Aligned::OccTable<Sigma>>{}, "compact2Aligned");
//    cb(std::type_identity<occtable::wavelet::OccTable<Sigma>>{}, "wavelet");
//    cb(std::type_identity<occtable::sdsl_wt_bldc::OccTable<Sigma>>{}, "sdsl_wt_bldc");
//    cb(std::type_identity<occtable::sdsl_wt_epr::OccTable<Sigma>>{}, "sdsl_wt_epr");
//    cb(std::type_identity<occtable::compact::OccTable<Sigma>>{}, "compact");
//    cb(std::type_identity<occtable::compactAligned::OccTable<Sigma>>{}, "compactAligned");
//    cb(std::type_identity<occtable::compactPrefix::OccTable<Sigma>>{}, "compactPrefix");
//    cb(std::type_identity<occtable::bitvectorPrefix::OccTable<Sigma>>{}, "bitvectorPrefix");
}


struct Query {
    std::string name;
    bool reverse;
    Query(std::string _name, bool _reverse)
        : name{std::move(_name)}
        , reverse{_reverse}
        {}
};
template <size_t Sigma>
auto loadQueries(std::string path, bool reverse) {
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
                if (*ptr != '>') {
                    throw std::runtime_error("expected '>'");
                }
                ++ptr;
                if (*ptr == ' ') {
                    ++ptr;
                }
                while (ptr != (b.data() + b.size()) and *ptr != '\n') {
                    name += *ptr;
                    ++ptr;
                }
                ++ptr;
                mode = Mode::Sequence;
                queryInfos.emplace_back(name, false);
                if (reverse) {
                    queryInfos.emplace_back(name, true);
                }
            } else if (mode == Mode::Sequence) {
                if (*ptr == '>' || (ptr+1) == (b.data() + b.size())) {
                    queries.push_back(query);
                    if (reverse) {
                        std::reverse(query.begin(), query.end());
                        for (auto& c : query) {
                            if (c == 1) c = 4;
                            else if (c == 2) c = 3;
                            else if (c == 3) c = 2;
                            else if (c == 4) c = 1;
                        }
                        queries.push_back(query);
                    }
                    query.clear();
                    mode = Mode::Name;
                    if ((ptr+1) == (b.data() + b.size())) {
                        ++ptr;
                    }
                } else {
                    if (*ptr == '$') query.push_back(0);
                    else if (*ptr == 'A' || *ptr == 'a') query.push_back(1);
                    else if (*ptr == 'C' || *ptr == 'c') query.push_back(2);
                    else if (*ptr == 'G' || *ptr == 'g') query.push_back(3);
                    else if (*ptr == 'T' || *ptr == 't') query.push_back(4);
                    else if ((*ptr == 'N' || *ptr == 'n') and Sigma == 6) query.push_back(5);
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

template <size_t Sigma, typename CSA, template <size_t> typename Table>
auto loadIndex(std::string path) {
    auto indexPath = path + "." + Table<Sigma>::extension() + ".index";
    if (!std::filesystem::exists(indexPath)) {
        auto [ref, refInfo] = loadQueries<Sigma>(path + ".fasta", false);
        auto index = BiFMIndex<Table<Sigma>>{ref, 16};
        // save index here
        auto ofs     = std::ofstream{indexPath, std::ios::binary};
        auto archive = cereal::BinaryOutputArchive{ofs};
        archive(index);
        return index;
    } else {
        auto ifs     = std::ifstream{indexPath, std::ios::binary};
        auto archive = cereal::BinaryInputArchive{ifs};
        auto index = BiFMIndex<Table<Sigma>>{cereal_tag{}};
        archive(index);
        return index;
    }
}




int main(int argc, char const* const* argv) {
    StopWatch watch;
    constexpr size_t Sigma = 5;


//    std::string algorithm = "pseudo";
    std::string generator = "h2";
    size_t maxQueries{};
    size_t readLength{};
    bool saveOutput{false};
    size_t minK{0}, maxK{6}, k_stepSize{1};
    bool reverse{true};

    std::vector<std::string> algorithms;

    for (int i{1}; i < argc; ++i) {
        if (argv[i] == std::string{"--algo"} and i+1 < argc) {
            ++i;
            algorithms.emplace_back(argv[i]);
//            algorithm = argv[i];
        } else if (argv[i] == std::string{"--gen"} and i+1 < argc) {
            ++i;
            generator = argv[i];
        } else if (argv[i] == std::string{"--queries"} and i+1 < argc) {
            ++i;
            maxQueries = std::stod(argv[i]);
        } else if (argv[i] == std::string{"--read_length"} and i+1 < argc) {
            ++i;
            readLength = std::stod(argv[i]);
        } else if (argv[i] == std::string{"--save_output"}) {
            saveOutput = true;
        } else if (argv[i] == std::string{"--min_k"} and i+1 < argc) {
            ++i;
            minK = std::stod(argv[i]);
        } else if (argv[i] == std::string{"--max_k"} and i+1 < argc) {
            ++i;
            maxK = std::stod(argv[i]);
        } else if (argv[i] == std::string{"--stepSize_k"} and i+1 < argc) {
            ++i;
            k_stepSize = std::stod(argv[i]);
        } else if (argv[i] == std::string{"--no-reverse"}) {
            reverse = false;
        } else {
            throw std::runtime_error("unknown commandline " + std::string{argv[i]});
        }
    }
     auto const [queries, queryInfos] = loadQueries<Sigma>("/home/gene/hg38/sampled_illumina_reads.fasta", reverse);

//    auto const [queries, queryInfos] = loadQueries<Sigma>("/home/gene/short_test/read.fasta", reverse);

    fmt::print("loaded {} queries (incl reverse complements)\n", queries.size());

    fmt::print("{:15}: {:>10}  ({:>10} +{:>10} ) {:>10}    - results: {:>10}/{:>10}/{:>10}/{:>10} - mem: {:>13}\n", "name", "time_search + time_locate", "time_search", "time_locate", "(time_search+time_locate)/queries.size()", "resultCt", "results.size()", "uniqueResults.size()", "readIds.size()", "memory");



//    auto [bwt, bwtRev] = generateBWT<Sigma>(1ul<<20);
    visitAllTables<Sigma>([&]<template <size_t> typename Table>(std::type_identity<Table<Sigma>>, std::string name) {
        if constexpr (OccTableMemoryUsage<Table<Sigma>>) {
            size_t s = Table<Sigma>::expectedMemoryUsage(3'000'000'000ul);
    //        fmt::print("expected memory: {}\n", s);
            if (s > 1024*1024*1024*56ul) {
                fmt::print("{} skipping, to much memory\n", name);
                return;
            }
        }

        auto index = loadIndex<Sigma, CSA, Table>("/home/gene/hg38/text.dna4");

//        auto index = loadIndex<Sigma, CSA, Table>("/home/gene/short_test/text");
        fmt::print("index loaded\n");
        for (auto const& algorithm : algorithms) {
            fmt::print("using algorithm {}\n", algorithm);

            auto mut_queries = queries;
            if (maxQueries != 0) {
                mut_queries.resize(std::min(mut_queries.size(), maxQueries));
            }
            if (readLength != 0) {
                for (auto& q : mut_queries) {
                    q.resize(std::min(readLength, q.size()));
                }
            }



            auto memory = [&] () -> size_t {
                if constexpr (OccTableMemoryUsage<Table<Sigma>>) {
                    return index.memoryUsage();
                } else {
                    return 0ul;
                }
            }();
            for (size_t k{minK}; k<=maxK; k = k + k_stepSize) {
                if (k >= 4 and k != 6) {
                    mut_queries.resize(mut_queries.size() / 10);
                }
                auto search_scheme = [&]() {
                    if (generator == "pigeon_opt")   return search_schemes::expand(search_schemes::generator::pigeon_opt(0, k), mut_queries[0].size());
                    if (generator == "pigeon")       return search_schemes::expand(search_schemes::generator::pigeon_trivial(0, k), mut_queries[0].size());
                    if (generator == "h2")           return search_schemes::expand(search_schemes::generator::h2(k+2, 0, k), mut_queries[0].size());
                    if (generator == "kucherov")     return search_schemes::expand(search_schemes::generator::kucherov(k+1, k), mut_queries[0].size());
                    if (generator == "backtracking") return search_schemes::expand(search_schemes::generator::backtracking(1, 0, k), mut_queries[0].size());
                    throw std::runtime_error("unknown search scheme");
                }();

    //            auto search_scheme = search_schemes::expand(search_schemes::generator::pigeon_opt(0, k), mut_queries[0].size());
    //            auto search_scheme = search_schemes::expand(search_schemes::generator::pigeon_trivial(0, k), mut_queries[0].size());
    //            auto search_scheme = search_schemes::expand(search_schemes::generator::h2(k+2, 0, k), mut_queries[0].size());
    //            auto search_scheme = search_schemes::expand(search_schemes::generator::kucherov(k+1, 0, k), mut_queries[0].size());
    //            auto search_scheme = search_schemes::expand(search_schemes::generator::backtracking(1, 0, k), mut_queries[0].size());

        //
                size_t resultCt{};
                StopWatch sw;
                std::vector<std::tuple<size_t, size_t, size_t, std::string>> results{};
                std::vector<std::tuple<size_t, LeftBiFMIndexCursor<decltype(index)>, size_t, std::string>> resultCursors;

                auto res_cb = [&](size_t queryId, auto cursor, size_t errors) {
                    resultCursors.emplace_back(queryId, cursor, errors, "");
                };
                auto res_cb2 = [&](size_t queryId, auto cursor, size_t errors, auto const& actions) {
                    std::string s;
                    for (auto a : actions) {
                        s += a;
                    }
                    resultCursors.emplace_back(queryId, cursor, errors, s);
                };



                if (algorithm == "pseudo") search_pseudo::search<true>(index, mut_queries, search_scheme, res_cb);
                else if (algorithm.size() == 15 && algorithm.substr(0, 13) == "pseudo_fmtree")  search_pseudo::search<true>(index, mut_queries, search_scheme, res_cb);
                else if (algorithm == "pseudo_fmtree")  search_pseudo::search<true>(index, mut_queries, search_scheme, res_cb);
                else if (algorithm == "ng12") search_ng12::search(index, mut_queries, search_scheme, res_cb);
                else if (algorithm == "ng14") search_ng14::search(index, mut_queries, search_scheme, res_cb);
                else if (algorithm == "ng15") search_ng15::search(index, mut_queries, search_scheme, res_cb);
                else if (algorithm == "ng16") search_ng16::search(index, mut_queries, search_scheme, res_cb);
                else if (algorithm == "ng17") search_ng17::search(index, mut_queries, search_scheme, res_cb);
                else if (algorithm == "ng20") search_ng20::search(index, mut_queries, search_scheme, res_cb);
                else if (algorithm == "ng21") search_ng21::search(index, mut_queries, search_scheme, res_cb);
                else if (algorithm == "ng22") search_ng22::search(index, mut_queries, search_scheme, res_cb2);
                else if (algorithm == "noerror") search_no_errors::search(index, mut_queries, res_cb);

    //            search_ng14::search(index, mut_queries, search_scheme, [&](size_t queryId, auto cursor, size_t errors) {
    //            search_ng15::search(index, mut_queries, search_scheme, [&](size_t queryId, auto const& cursor, size_t errors) {
    //            search_ng16::search(index, mut_queries, search_scheme, [&](size_t queryId, auto const& cursor, size_t errors) {
    //            search_ng17::search(index, mut_queries, search_scheme, [&](size_t queryId, auto const& cursor, size_t errors) {
    //            search_ng20::search(index, mut_queries, search_scheme, [&](size_t queryId, auto cursor, size_t errors) {
    //            search_ng21::search(index, mut_queries, search_scheme, [&](size_t queryId, auto cursor, size_t errors) {

    //                resultCursors.emplace_back(queryId, cursor, errors);
    //            });

                auto time_search = sw.reset();

                if (algorithm.size() == 15 && algorithm.substr(0, 13)  == "pseudo_fmtree") {
                    size_t maxDepth = std::stod(algorithm.substr(13, 2));
                    for (auto const& [queryId, cursor, e, action] : resultCursors) {
                        for (auto [seqId, pos] : LocateFMTree{index, cursor, maxDepth}) {
                            results.emplace_back(queryId, pos, e, action);
                        }
                        resultCt += cursor.len;
                    }
                } else if (algorithm == "pseudo_fmtree") {
                    for (auto const& [queryId, cursor, e, action] : resultCursors) {
                        locateFMTree<16>(index, cursor, [&](size_t seqId, size_t pos) {
                            (void)seqId;
                            results.emplace_back(queryId, pos, e ,action);
                        });
                        resultCt += cursor.len;
                    }

                } else {
                    for (auto const& [queryId, cursor, e, action] : resultCursors) {
                        for (auto [seqId, pos] : LocateLinear{index, cursor}) {
                            results.emplace_back(queryId, pos, e, action);
                        }
                        resultCt += cursor.len;
                    }
                }
                auto time_locate = sw.reset();

                auto uniqueResults = [](auto list) {
                    std::sort(begin(list), end(list));
                    list.erase(std::unique(begin(list), end(list)), list.end());
                    return list;
                }(results);
                std::unordered_set<size_t> readIds;
                for (auto const& [queryId, cursor, e, actions] : resultCursors) {
                    if (queryId > mut_queries.size()/2) {
                        readIds.insert(queryId - mut_queries.size() / 2);
                    } else {
                        readIds.insert(queryId);
                    }
                }

                fmt::print("{:15} {:3}: {:>10.3}s ({:>10.3}s+{:>10.3}s) {:>10.3}q/s - results: {:>10}/{:>10}/{:>10}/{:>10} - mem: {:>13}\n", name, k, time_search + time_locate, time_search, time_locate, mut_queries.size() / (time_search+time_locate), resultCt, results.size(), uniqueResults.size(), readIds.size(), memory);
                {
                    if (saveOutput) {
                        auto filename =fmt::format("out.k{}.alg{}.ss{}.txt", k, algorithm, name);
                        auto ofs = fopen(filename.c_str(), "w");
                        fmt::print(ofs, "identifier\tposition\tlength\tED\treverseComlpement\talignment\n");
                        for (auto const& [queryId, pos, e, action] : results) {
                            auto const& qi = queryInfos[queryId];
                            fmt::print(ofs, "{}\t{}\t{}\t{}\t{}\t{}\n", qi.name, pos, 101, e, qi.reverse?1:0, action);
                        }
                        fclose(ofs);
                    }
                }
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


