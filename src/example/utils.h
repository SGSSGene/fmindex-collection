// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "utils/utils.h"

#include <fmindex-collection/fmindex-collection.h>
#include <fmindex-collection/fmindex/merge.h>
#include <fmindex-collection/suffixarray/DenseCSA.h>

#include <cereal/archives/binary.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/vector.hpp>
#include <string>
#include <vector>

struct Query {
    std::string name;
    bool reverse;
    Query(std::string _name, bool _reverse)
        : name{std::move(_name)}
        , reverse{_reverse}
        {}
};
template <size_t Sigma>
auto loadQueries(std::string path, bool reverse, bool convertUnknownChar) {
    std::vector<std::vector<uint8_t>> queries;
    std::vector<Query> queryInfos;
    if (path.empty() || !std::filesystem::exists(path)) {
        return std::make_tuple(queries, queryInfos);
    }
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
                        std::ranges::reverse(query);
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
                        if (convertUnknownChar) {
                            if (Sigma == 6) {
                                query.push_back(5);
                            } else {
                                query.push_back(1);
                            }
                        } else {
                            throw std::runtime_error("unknown alphabet");
                        }
                    }
                    ++ptr;
                }
            }
        }
    }
    return std::make_tuple(queries, queryInfos);
}

template <typename CSA, size_t Sigma, template <size_t> typename String>
auto loadIndex(std::string path, size_t samplingRate, size_t threadNbr, bool convertUnknownChar) {
    auto sw = StopWatch{};
    auto indexPath = path + ".index";
    if (!std::filesystem::exists(indexPath)) {
        auto [ref, refInfo] = loadQueries<Sigma>(path, false, convertUnknownChar);
        auto index = [&]() {
            auto refs = std::vector<std::vector<uint8_t>>{};
            refs.resize(1);

            auto index = std::optional<fmc::BiFMIndex<Sigma, String>>{};
            for (auto& r : ref) {
                refs[0] = std::move(r);
                auto newIndex = fmc::BiFMIndex<Sigma, String>{refs, samplingRate, threadNbr};
                if (!index) {
                    index = std::move(newIndex);
                } else {
                    index = fmc::fmindex::merge(*index, newIndex);
                }
            }
            return *index;
        }();
        // save index here
        auto ofs     = std::ofstream{indexPath, std::ios::binary};
        auto archive = cereal::BinaryOutputArchive{ofs};
        archive(index);
        return index;
    } else {
        auto ifs     = std::ifstream{indexPath, std::ios::binary};
        auto archive = cereal::BinaryInputArchive{ifs};
        auto index = fmc::BiFMIndex<Sigma, String>{};
        archive(index);
        std::cout << "loading took " << sw.peek() << "s\n";
        return index;
    }
}

template <typename CSA, size_t Sigma, template <size_t> typename String>
auto loadDenseIndex(std::string path, size_t samplingRate, size_t threadNbr, bool partialBuildUp, bool convertUnknownChar) {
    auto sw = StopWatch{};
    auto indexPath = path + ".tab.dense.index";

    if (!std::filesystem::exists(indexPath)) {
        auto [ref, refInfo] = loadQueries<Sigma>(path, false, convertUnknownChar);
        using Index = fmc::BiFMIndex<Sigma, String>;
        auto index = [&]() -> Index {
            if (!partialBuildUp) {
                return {ref, samplingRate, threadNbr};
            }
            auto longestRef = std::accumulate(ref.begin(), ref.end(), size_t{}, [](size_t a, auto const& v) {
                return std::max(a, v.size());
            });
            std::cout << "longest ref: " << longestRef << "\n";
            #if 0
            auto refs = std::vector<std::vector<uint8_t>>{};
            refs.resize(1);
            auto index = std::optional<Index>{};
            for (auto& r : ref) {
                refs[0] = std::move(r);
                std::cout << "indexing " << refs[0].size() << "\n";
                auto newIndex = Index{refs, samplingRate, threadNbr};
                if (!index) {
                    index = std::move(newIndex);
                } else {
                    std::cout << "merging " << index->size() << " + " << refs[0].size() << "\n";
                    index = merge(*index, newIndex);
                }
            }
            return std::move(*index);
            #elif 1
            auto refs = std::vector<std::vector<uint8_t>>{};
            auto index = std::optional<Index>{};
            size_t acc = 0;
            auto makePartialIndex = [&]() {
                if (refs.empty()) return;
                std::cout << "indexing " << acc << "\n";
                auto newIndex = Index{refs, samplingRate, threadNbr};
                if (!index) {
                    index = std::move(newIndex);
                } else {
                    std::cout << "merging " << index->size() << " + " << acc << "\n";
                    index = fmc::fmindex::merge(*index, newIndex);
                }

                acc = 0;
                refs.clear();
            };
            for (size_t i{}; i < ref.size(); ++i) {
                if (ref[i].size() + acc >= longestRef) {
                    makePartialIndex();
                }
                refs.emplace_back(std::move(ref[i]));
                acc += refs.back().size();
            }
            makePartialIndex();
            return std::move(*index);

            #else
            auto refs = std::vector<std::vector<uint8_t>>{};
            refs.resize(1);
            auto indices = std::vector<Index>{};

            auto sort = [&]() {
                std::sort(indices.begin(), indices.end(), [](auto const& lhs, auto const& rhs) {
                    return lhs.size() > rhs.size();
                });
            };
            for (auto& r : ref) {
                refs[0] = std::move(r);
                std::cout << "indexing " << refs[0].size() << " " << indices.size() << "\n";
                auto newIndex = Index{refs, samplingRate, threadNbr};
                indices.emplace_back(std::move(newIndex));
                sort();

                while (indices.size() > 1) {
                    auto const& l1 = *(indices.end()-1);
                    auto const& l2 = *(indices.end()-2);
                    if (l2.size() > l1.size()*2) {
                        break;
                    }
                    std::cout << "merging " << l2.size() << " + " << l1.size() << " " << indices.size() << "\n";
                    auto newIndex = fmc::fmindex::merge(l2, l1);
                    indices.pop_back(); indices.pop_back();
                    indices.emplace_back(std::move(newIndex));
                    sort();
                }
            }
            while (indices.size() > 1) {
                auto const& l1 = *(indices.end()-1);
                auto const& l2 = *(indices.end()-2);
                std::cout << "merging " << l2.size() << " + " << l1.size() << " " << indices.size() << "(fin)\n";
                auto newIndex = fmc::fmindex::merge(l2, l1);
                indices.pop_back(); indices.pop_back();
                indices.emplace_back(std::move(newIndex));
                sort();
            }
            return std::move(indices.back());
            #endif
        }();

        // save index here
        auto ofs     = std::ofstream{indexPath, std::ios::binary};
        auto archive = cereal::BinaryOutputArchive{ofs};
        archive(index);
        return index;
    } else {
        auto ifs     = std::ifstream{indexPath, std::ios::binary};
        auto archive = cereal::BinaryInputArchive{ifs};
        auto index = fmc::BiFMIndex<Sigma, String>{};
        archive(index);
        std::cout << "loading took " << sw.peek() << "s\n";
        return index;
    }
}

template <size_t Sigma, typename CB>
void visitAllStrings(CB cb) {
    cb.template operator()<Sigma, fmc::string::InterleavedBitvector16>();
}
