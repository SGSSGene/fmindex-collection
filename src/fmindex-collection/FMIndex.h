#pragma once

#include "occtable/concepts.h"

namespace fmindex_collection {

template <OccTable Table>
struct FMIndex {
    static size_t constexpr Sigma = Table::Sigma;

    Table occTable;

    FMIndex(std::vector<uint8_t> const& _bwt)
        : occTable{_bwt}
    {}

    size_t size() const {
        return occTable.size();
    }
    size_t rank(size_t idx, uint8_t symb) const {
        return occTable.rank(idx, symb);
    }
    size_t prefix_rank(size_t idx, uint8_t symb) const {
        return occTable.prefix_rank(idx, symb);
    }
};

}
