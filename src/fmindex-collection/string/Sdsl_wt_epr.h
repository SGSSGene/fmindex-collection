// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#ifdef FMC_USE_SDSL

#include "concepts.h"
#include "utils.h"

#include <array>
#include <cstdint>
#include <filesystem>
#include <sdsl/cereal.hpp>
#include <sdsl/construct.hpp>
#include <sdsl/suffix_trees.hpp>
#include <sdsl/wt_epr.hpp>
#include <span>
#include <vector>

namespace fmc::string {
/** !TODO This structure is not working as expected
 */
template <size_t TSigma>
struct Sdsl_wt_epr {
    static constexpr size_t Sigma = TSigma;

    using sdsl_wt_index_type = sdsl::wt_epr<TSigma>;
    sdsl_wt_index_type index;
    size_t             totalLength{};

    Sdsl_wt_epr() = default;
    Sdsl_wt_epr(std::span<uint8_t const> in_symbols) {
        totalLength = in_symbols.size();
        auto _symbols = std::vector<uint8_t>{begin(in_symbols), end(in_symbols)};
        {
            auto ofs = std::ofstream{"tmp.sdsl.tmp", std::ios::binary};
            ofs.write(reinterpret_cast<char const*>(_symbols.data()), _symbols.size());
        }

        sdsl::construct(index, "tmp.sdsl.tmp", 1);
    }

    size_t size() const {
        return totalLength;
    }

    uint8_t symbol(uint64_t idx) const {
        assert(idx < totalLength);

        assert([&]() {
            for (size_t s{0}; s < Sigma; ++s) {
                auto diff = rank(idx+1, s) - rank(idx, s);
                if (diff != 0) {
                    return index[idx] == s;
                }
            }
            return true;
        }());

        return index[idx];
    }

    uint64_t rank(uint64_t idx, uint8_t symb) const {
        assert(idx <= totalLength);
        assert(symb < TSigma);

        return index.rank(idx, symb);
    }

    uint64_t prefix_rank(uint64_t idx, uint8_t symb) const {
        assert(idx <= totalLength);
        assert(symb <= TSigma);

        uint64_t a{};
        for (size_t i{0}; i < symb; ++i) {
            a += index.rank(idx, i);
        }
        return a;
    }

    auto all_ranks(uint64_t idx) const -> std::array<uint64_t, Sigma> {
        assert(idx <= totalLength);

        std::array<uint64_t, Sigma> rs{0};
        for (size_t i{0}; i < Sigma; ++i) {
            rs[i] = rank(idx, i);
        }
        return rs;
    }

    auto all_ranks_and_prefix_ranks(uint64_t idx) const -> std::tuple<std::array<uint64_t, Sigma>, std::array<uint64_t, Sigma>> {
        assert(idx <= totalLength);

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
        ar(index, totalLength);
    }
};
//(Sdsl wt epr seems broken?, since it doesn't work for Sigma <= 2 and Sigma == 256
//static_assert(checkString_c<Sdsl_wt_epr>);

}
#endif
