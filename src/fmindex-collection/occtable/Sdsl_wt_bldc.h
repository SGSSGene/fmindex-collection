// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include "concepts.h"

#include <array>
#include <cstdint>
#include <filesystem>
#include <sdsl/construct.hpp>
#include <sdsl/suffix_trees.hpp>
#include <span>
#include <vector>

namespace fmindex_collection {
namespace occtable {
namespace sdsl_wt_bldc {

inline void writeFile(std::filesystem::path const& file, std::vector<uint8_t> const& buffer) {
    auto ofs = std::ofstream{file, std::ios::binary};
    ofs.write(reinterpret_cast<char const*>(buffer.data()), buffer.size());
}

using sdsl_wt_index_type =
    sdsl::wt_blcd<sdsl::bit_vector, // Wavelet tree type
                    sdsl::rank_support_v<1>,
                    sdsl::select_support_scan<>,
                    sdsl::select_support_scan<0>>;

template <size_t TSigma>
struct OccTable {
    using TLengthType = uint64_t;
    static constexpr size_t Sigma = TSigma;

    sdsl_wt_index_type index;


    std::array<uint64_t, Sigma+1> C{};

    OccTable(std::span<uint8_t const> in_bwt) {
        auto _bwt = std::vector<uint8_t>{begin(in_bwt), end(in_bwt)};
        for (auto& c : _bwt) {
            c += 1;
        }
        writeFile("tmp.sdsl.tmp", _bwt);
        sdsl::construct(index, "tmp.sdsl.tmp", 1);

        C[0] = 0;
        for (size_t i{0}; i < Sigma; ++i) {
            auto r = index.rank(_bwt.size(), i+1);
            C[i+1] = r + C[i];
        }
    }

    OccTable(cereal_tag) {}

    static auto name() -> std::string {
        return "SDSL Wavelet";
    }

    static auto extension() -> std::string {
        return "sdslwt";
    }

    uint64_t size() const {
        return C.back();
    }

    uint64_t rank(uint64_t idx, uint8_t symb) const {
        return index.rank(idx, symb+1) + C[symb];
    }

    uint64_t prefix_rank(uint64_t idx, uint8_t symb) const {
        uint64_t a{};
        for (size_t i{0}; i <= symb; ++i) {
            a += index.rank(idx, i+1);
        }
        return a;
    }

    size_t symbol(uint64_t idx) const {
        idx += 1;
        for (size_t i{0}; i < Sigma-1; ++i) {
            if (rank(idx, i) > rank(idx-1, i)) {
                return i;
            }
        }
        return Sigma-1;
    }

    auto all_ranks(uint64_t idx) const -> std::tuple<std::array<uint64_t, Sigma>, std::array<uint64_t, Sigma>> {
        std::array<uint64_t, Sigma> rs{0};
        std::array<uint64_t, Sigma> prs{0};
        for (size_t i{0}; i < Sigma; ++i) {
            rs[i] = rank(idx, i);
            prs[i] = prefix_rank(idx, i);
        }
        return {rs, prs};
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(index, C);
    }
};
static_assert(checkOccTable<OccTable>);

}
}
}
