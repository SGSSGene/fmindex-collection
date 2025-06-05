// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause

#include <cereal/archives/binary.hpp>
#include <filesystem>
#include <fstream>

namespace fmindex_collection {

// loads an fm index from disk
template <typename Index>
void saveIndex(Index const& _index, std::filesystem::path _fileName) {
    auto ofs     = std::ofstream(_fileName, std::ios::binary);
    auto archive = cereal::BinaryOutputArchive{ofs};
    archive(_index);
}

// loads an fm index from disk
template <typename Index>
auto loadIndex(std::filesystem::path _fileName) {
    auto ifs     = std::ifstream(_fileName, std::ios::binary);
    auto archive = cereal::BinaryInputArchive{ifs};
    auto index = Index{};
    archive(index);
    return index;
}

}
