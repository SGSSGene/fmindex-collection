#include "occtable/concepts.h"

template <OccTable Table>
struct BiFMIndex {
    static size_t constexpr Sigma = Table::Sigma;

    Table occ;
    Table occRev;

    BiFMIndex(std::vector<uint8_t> const& bwt, std::vector<uint8_t> const& bwtRev)
        : occ{bwt}
        , occRev{bwtRev}
    {
        assert(bwt.size() == bwtRev.size());
        assert(occ.size() == occRev.size());
    }

    size_t size() const {
        return occ.size();
    }
};

template <typename Index>
struct BiFMIndexCursor {
    static size_t constexpr Sigma = Index::Sigma;

    Index const* index;
    size_t lb;
    size_t lbRev;
    size_t len;
    BiFMIndexCursor(Index const& index)
        : BiFMIndexCursor{index, 0, 0, index.size()}
    {}
    BiFMIndexCursor(Index const& index, size_t lb, size_t lbRev, size_t len)
        : index{&index}
        , lb{lb}
        , lbRev{lbRev}
        , len{len}
    {}
    auto extendLeft() const -> std::array<BiFMIndexCursor, Sigma> {
        auto const& occ = index->occ;
        auto [rs1, prs1] = occ.all_ranks(lb);
        auto [rs2, prs2] = occ.all_ranks(lb+len);

        auto cursors = std::array<BiFMIndexCursor, Sigma>{*index};
        cursors[0] = BiFMIndexCursor{*index, rs1[0], lbRev, prs2[0] - rs1[0]};
        for (size_t i{1}; i < Sigma; ++i) {
            cursors[i] = BiFMIndexCursor{*index, rs1[i], lbRev + prs2[i-1] - prs1[i-1], prs2[i] - rs1[i]};
        }
        return cursors;
    }

    auto extendRight() const -> std::array<BiFMIndexCursor, Sigma> {
        auto const& occ = index->occ;
        auto [rs1, prs1] = occ.all_ranks(lbRev);
        auto [rs2, prs2] = occ.all_ranks(lbRev+len);

        auto cursors = std::array<BiFMIndexCursor, Sigma>{*index};
        cursors[0] = BiFMIndexCursor{*index, lb, rs1[0], prs2[0] - rs1[0]};
        for (size_t i{1}; i < Sigma; ++i) {
            cursors[i] = BiFMIndexCursor{*index, lbRev + prs2[i-1] - prs1[i-1], rs1[i], prs2[i] - rs1[i]};
        }
        return cursors;
    }

    auto extendLeft(uint8_t symb) const -> BiFMIndexCursor {
        assert(symb > 0);
        auto& occ = index->occ;
        size_t newLb    = occ.rank(lb, symb);
        size_t newLbRev = lbRev + occ.prefix_rank(lb+len, symb-1) - occ.prefix_rank(lb, symb-1);
        size_t newLen   = occ.rank(lb+len, symb) - newLb;
        return {*index, newLb, newLbRev, newLen};
    }
    auto extendRight(uint8_t symb) const -> BiFMIndexCursor {
        assert(symb > 0);
        auto& occ = index->occRev;
        size_t newLb    = lb + occ.prefix_rank(lbRev+len, symb-1) - occ.prefix_rank(lbRev, symb-1);
        size_t newLbRev = occ.rank(lbRev, symb);
        size_t newLen   = occ.rank(lbRev+len, symb) - newLbRev;

        return {*index, newLb, newLbRev, newLen};
    }

};
