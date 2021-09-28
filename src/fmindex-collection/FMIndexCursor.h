#pragma once

#include "FMIndex.h"

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
