#pragma once

#include "concepts.h"

#include <array>
#include <cstdint>
#include <vector>

namespace occtable {
namespace naive {

template <size_t TSigma>
struct OccTable {
    static constexpr size_t Sigma = TSigma;

    std::vector<std::array<uint64_t, Sigma>> occ{};
    std::array<uint64_t, Sigma+1> C{};

    static size_t expectedMemoryUsage(size_t length) {
        size_t C       = sizeof(uint64_t) * (Sigma+1);
        size_t entries = sizeof(uint64_t) * Sigma * (length+1);
        return C + entries;
    }

    OccTable(std::vector<uint8_t> const& _bwt) {
        occ.reserve(_bwt.size()+1);
        occ.push_back(std::array<uint64_t, Sigma>{});

        for (size_t i{0}; i < _bwt.size(); ++i) {
            occ.push_back(occ.back());
            occ.back()[_bwt[i]] += 1;
        }
        C[0] = 0;
        for (size_t i{0}; i < Sigma; ++i) {
            C[i+1] = occ.back()[i] + C[i];
        }
    }

    uint64_t size() const {
        return C.back();
    }

    uint64_t rank(uint8_t symb, uint64_t idx) const {
        return occ[idx][symb] + C[symb];
    }

    uint64_t prefix_rank(uint8_t symb, uint64_t idx) const {
        uint64_t a{};
        for (size_t i{0}; i <= symb; ++i) {
            a += occ[idx][i];
        }
        return a;
    }
};

static_assert(checkOccTable<OccTable>);

}
}
