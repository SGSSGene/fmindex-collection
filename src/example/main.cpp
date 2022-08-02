#include "utils.h"
#include "argp.h"

#include <fmindex-collection/fmindex-collection.h>
#include <fmindex-collection/occtable/all.h>
#include <fmindex-collection/search/all.h>
#include <fmindex-collection/locate.h>
#include <search_schemes/generator/all.h>
#include <search_schemes/expand.h>

#include <fmt/format.h>
#include <unordered_set>

using namespace fmindex_collection;

template <size_t Sigma, typename CB>
void visitAllTables(CB cb) {
//    cb(std::type_identity<occtable::naive::OccTable<Sigma>>{}, "naive");
//    cb(std::type_identity<occtable::bitvector::OccTable<Sigma>>{}, "bitvector");
//    cb(std::type_identity<occtable::interleaved16::OccTable<Sigma>>{}, "interleaved16");
//    cb(std::type_identity<occtable::interleaved32::OccTable<Sigma>>{}, "interleaved32");
    cb(std::type_identity<occtable::interleavedEPR16::OccTable<Sigma>>{});
    cb(std::type_identity<occtable::interleavedEPR32::OccTable<Sigma>>{});
    cb(std::type_identity<occtable::interleavedEPR16V2::OccTable<Sigma>>{});
    cb(std::type_identity<occtable::interleavedEPR32V2::OccTable<Sigma>>{});
    cb(std::type_identity<occtable::eprV4::OccTable<Sigma>>{});
    cb(std::type_identity<occtable::eprV5::OccTable<Sigma>>{});
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


int main(int argc, char const* const* argv) {
    constexpr size_t Sigma = 5;


    auto config = loadConfig(argc, argv);
    if (config.help) {
        fmt::print("Usage:\n"
                    "./example --index somefile.fasta\n"
                    "   this will only build the index files for somefile.fasta\n"
                    "\n"
                    "./example --index somefile.fasta\\\n"
                    "          --query queryfile.fasta\\\n"
                    "          --algo [pseudo, pseudo_fmtree00-pseudo_fmtree99, ng12, ng14, ng15, ng16, ng17, ng20, ng21, ng22, noerror]\\\n"
                    "          --gen <pigeon_opt|pigeon|h2|kucherov|backtracking>\\\n"
                    "          --queries <int> (maximal of number of queries)\\\n"
                    "          --read_length <int> (shorten all queries to this length)\\\n"
                    "          --save_output (saves output at the end)\\\n"
                    "          --min_k <int> (minimal number of errors)\\\n"
                    "          --max_k <int> (maximal number of errors)\\\n"
                    "          --stepSize_k <int> (steps of errors)\\\n"
                    "          --no-reverse (don't use reverse compliment)\n"
        );
        return 0;
    }
    auto const [queries, queryInfos] = loadQueries<Sigma>(config.queryPath, config.reverse);

    if (!queries.empty()) {
        fmt::print("loaded {} queries (incl reverse complements)\n", queries.size());
        fmt::print("{:15}: {:>10}  ({:>10} +{:>10} ) {:>10}    - results: {:>10}/{:>10}/{:>10}/{:>10} - mem: {:>13}\n", "name", "time_search + time_locate", "time_search", "time_locate", "(time_search+time_locate)/queries.size()", "resultCt", "results.size()", "uniqueResults.size()", "readIds.size()", "memory");
    }


    visitAllTables<Sigma>([&]<template <size_t> typename Table>(std::type_identity<Table<Sigma>>) {
        std::string name = Table<Sigma>::extension();

        if constexpr (OccTableMemoryUsage<Table<Sigma>>) {
            size_t s = Table<Sigma>::expectedMemoryUsage(3'000'000'000ul);
            if (s > 1024*1024*1024*56ul) {
                fmt::print("{} skipping, to much memory\n", name);
                return;
            }
        }

        auto index = loadIndex<Sigma, CSA, Table>(config.indexPath);
        fmt::print("index loaded\n");
        for (auto const& algorithm : config.algorithms) {
            fmt::print("using algorithm {}\n", algorithm);

            auto mut_queries = queries;
            if (config.maxQueries != 0) {
                mut_queries.resize(std::min(mut_queries.size(), config.maxQueries));
            }
            if (config.readLength != 0) {
                for (auto& q : mut_queries) {
                    q.resize(std::min(config.readLength, q.size()));
                }
            }



            auto memory = [&] () -> size_t {
                if constexpr (OccTableMemoryUsage<Table<Sigma>>) {
                    return index.memoryUsage();
                } else {
                    return 0ul;
                }
            }();
            for (size_t k{config.minK}; k <= config.maxK; k = k + config.k_stepSize) {
                if (k >= 4 and k != 6) {
                    mut_queries.resize(mut_queries.size() / 10);
                }
                auto search_scheme = [&]() {
                    if (config.generator == "pigeon_opt")   return search_schemes::expand(search_schemes::generator::pigeon_opt(0, k), mut_queries[0].size());
                    if (config.generator == "pigeon")       return search_schemes::expand(search_schemes::generator::pigeon_trivial(0, k), mut_queries[0].size());
                    if (config.generator == "h2")           return search_schemes::expand(search_schemes::generator::h2(k+2, 0, k), mut_queries[0].size());
                    if (config.generator == "kucherov")     return search_schemes::expand(search_schemes::generator::kucherov(k+1, k), mut_queries[0].size());
                    if (config.generator == "backtracking") return search_schemes::expand(search_schemes::generator::backtracking(1, 0, k), mut_queries[0].size());
                    throw std::runtime_error("unknown search scheme");
                }();

                size_t resultCt{};
                StopWatch sw;
                auto results       = std::vector<std::tuple<size_t, size_t, size_t, std::string>>{};
                auto resultCursors = std::vector<std::tuple<size_t, LeftBiFMIndexCursor<decltype(index)>, size_t, std::string>>{};

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
                    if (config.saveOutput) {
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
}
