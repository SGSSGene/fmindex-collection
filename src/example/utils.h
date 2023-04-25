#pragma once

#include "utils/utils.h"

#include <fmindex-collection/fmindex-collection.h>
#include <fmindex-collection/occtable/all.h>
#include <fmindex-collection/DenseCSA.h>

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
                        throw std::runtime_error("unknown alphabet");
                    }
                    ++ptr;
                }
            }
        }
    }
    return std::make_tuple(queries, queryInfos);
}

template <typename CSA, typename Table>
auto loadIndex(std::string path, size_t samplingRate, size_t threadNbr) {
    auto sw = StopWatch{};
    auto indexPath = path + "." + Table::extension() + ".index";
    if (!std::filesystem::exists(indexPath)) {
        auto [ref, refInfo] = loadQueries<Table::Sigma>(path, false);
        auto index = fmindex_collection::BiFMIndex<Table>{ref, samplingRate, threadNbr};
        // save index here
        auto ofs     = std::ofstream{indexPath, std::ios::binary};
        auto archive = cereal::BinaryOutputArchive{ofs};
        archive(index);
        return index;
    } else {
        auto ifs     = std::ifstream{indexPath, std::ios::binary};
        auto archive = cereal::BinaryInputArchive{ifs};
        auto index = fmindex_collection::BiFMIndex<Table>{fmindex_collection::cereal_tag{}};
        archive(index);
        std::cout << "loading took " << sw.peek() << "s\n";
        return index;
    }
}

template <typename CSA, typename Table>
auto loadDenseIndex(std::string path, size_t samplingRate, size_t threadNbr) {
    auto sw = StopWatch{};
    auto indexPath = path + "." + Table::extension() + ".dense.index";
    if (!std::filesystem::exists(indexPath)) {
        auto [ref, refInfo] = loadQueries<Table::Sigma>(path, false);
        auto index = fmindex_collection::BiFMIndex<Table, fmindex_collection::DenseCSA>{ref, samplingRate, threadNbr};
        // save index here
        auto ofs     = std::ofstream{indexPath, std::ios::binary};
        auto archive = cereal::BinaryOutputArchive{ofs};
        archive(index);
        return index;
    } else {
        auto ifs     = std::ifstream{indexPath, std::ios::binary};
        auto archive = cereal::BinaryInputArchive{ifs};
        auto index = fmindex_collection::BiFMIndex<Table, fmindex_collection::DenseCSA>{fmindex_collection::cereal_tag{}};
        archive(index);
        std::cout << "loading took " << sw.peek() << "s\n";
        return index;
    }
}


template <size_t Sigma, typename CB>
void visitAllTables(CB cb) {
    cb.template operator()<fmindex_collection::occtable::naive::OccTable<Sigma>>();
    cb.template operator()<fmindex_collection::occtable::bitvector::OccTable<Sigma>>();
    cb.template operator()<fmindex_collection::occtable::compactBitvector::OccTable<Sigma>>();
    cb.template operator()<fmindex_collection::occtable::compactBitvectorPrefix::OccTable<Sigma>>();
    cb.template operator()<fmindex_collection::occtable::interleaved8::OccTable<Sigma>>();
    cb.template operator()<fmindex_collection::occtable::interleaved16::OccTable<Sigma>>();
    cb.template operator()<fmindex_collection::occtable::interleaved32::OccTable<Sigma>>();
    cb.template operator()<fmindex_collection::occtable::interleaved8Aligned::OccTable<Sigma>>();
    cb.template operator()<fmindex_collection::occtable::interleaved16Aligned::OccTable<Sigma>>();
    cb.template operator()<fmindex_collection::occtable::interleaved32Aligned::OccTable<Sigma>>();
    cb.template operator()<fmindex_collection::occtable::wavelet::OccTable<Sigma>>();
    cb.template operator()<fmindex_collection::occtable::interleavedWavelet::OccTable<Sigma>>();
    cb.template operator()<fmindex_collection::occtable::interleavedPrefix::OccTable<Sigma>>();
#ifdef FMC_USE_SDSL
    cb.template operator()<fmindex_collection::occtable::sdsl_wt_bldc::OccTable<Sigma>>();
    cb.template operator()<fmindex_collection::occtable::sdsl_wt_epr::OccTable<Sigma>>();
#endif
    cb.template operator()<fmindex_collection::occtable::interleavedEPR8::OccTable<Sigma>>();
    cb.template operator()<fmindex_collection::occtable::interleavedEPR16::OccTable<Sigma>>();
    cb.template operator()<fmindex_collection::occtable::interleavedEPR32::OccTable<Sigma>>();
    cb.template operator()<fmindex_collection::occtable::interleavedEPR8Aligned::OccTable<Sigma>>();
    cb.template operator()<fmindex_collection::occtable::interleavedEPR16Aligned::OccTable<Sigma>>();
    cb.template operator()<fmindex_collection::occtable::interleavedEPR32Aligned::OccTable<Sigma>>();
    cb.template operator()<fmindex_collection::occtable::interleavedEPR8V2::OccTable<Sigma>>();
    cb.template operator()<fmindex_collection::occtable::interleavedEPR16V2::OccTable<Sigma>>();
    cb.template operator()<fmindex_collection::occtable::interleavedEPR32V2::OccTable<Sigma>>();
    cb.template operator()<fmindex_collection::occtable::interleavedEPR8V2Aligned::OccTable<Sigma>>();
    cb.template operator()<fmindex_collection::occtable::interleavedEPR16V2Aligned::OccTable<Sigma>>();
    cb.template operator()<fmindex_collection::occtable::interleavedEPR32V2Aligned::OccTable<Sigma>>();
    cb.template operator()<fmindex_collection::occtable::epr8V3::OccTable<Sigma>>();
    cb.template operator()<fmindex_collection::occtable::epr16V3::OccTable<Sigma>>();
    cb.template operator()<fmindex_collection::occtable::epr32V3::OccTable<Sigma>>();
    cb.template operator()<fmindex_collection::occtable::eprV4::OccTable<Sigma>>();
    cb.template operator()<fmindex_collection::occtable::eprV5::OccTable<Sigma>>();
    cb.template operator()<fmindex_collection::occtable::eprV6::OccTable<Sigma>>();
    cb.template operator()<fmindex_collection::occtable::interleavedEPRV7::OccTable<Sigma>>();
    cb.template operator()<fmindex_collection::occtable::interleavedEPRV7b::OccTable<Sigma>>();
    cb.template operator()<fmindex_collection::occtable::eprV8::OccTable<Sigma>>();
}




