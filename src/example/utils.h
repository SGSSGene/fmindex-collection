#pragma once

#include "utils/utils.h"

#include <fmindex-collection/fmindex-collection.h>

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
        auto index = fmindex_collection::BiFMIndex<Table<Sigma>>{ref, 16};
        // save index here
        auto ofs     = std::ofstream{indexPath, std::ios::binary};
        auto archive = cereal::BinaryOutputArchive{ofs};
        archive(index);
        return index;
    } else {
        auto ifs     = std::ifstream{indexPath, std::ios::binary};
        auto archive = cereal::BinaryInputArchive{ifs};
        auto index = fmindex_collection::BiFMIndex<Table<Sigma>>{fmindex_collection::cereal_tag{}};
        archive(index);
        return index;
    }
}





