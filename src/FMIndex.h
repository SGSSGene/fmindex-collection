#include "occtable/concepts.h"

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

template <typename Index>
struct FMIndexCursor {
    Index const* index;
    size_t lb;
    size_t len;
    FMIndexCursor(Index const& index)
        : FMIndexCursor{index, 0, index.size()}
    {}
    FMIndexCursor(Index const& index, size_t lb, size_t len)
        : index{&index}
        , lb{lb}
        , len{len}
    {}
    FMIndexCursor(FMIndexCursor const&) = default;
    auto operator=(FMIndexCursor const&) -> FMIndexCursor& = default;

    auto extendLeft(uint8_t symb) const -> FMIndexCursor {
        size_t newLb  = index->rank(lb, symb);
        size_t newLen = index->rank(lb+len, symb) - newLb;
        return {*index, newLb, newLen};
    }
};
