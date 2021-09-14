#include "occtable/concepts.h"

#include "CSA.h"

template <OccTable Table>
struct BiFMIndex {
    static size_t constexpr Sigma = Table::Sigma;

    Table occ;
    Table occRev;
    CSA   csa;



    BiFMIndex(std::vector<uint8_t> const& bwt, std::vector<uint8_t> const& bwtRev, CSA _csa)
        : occ{bwt}
        , occRev{bwtRev}
        , csa{std::move(_csa)}
    {
        assert(bwt.size() == bwtRev.size());
        assert(occ.size() == occRev.size());
        if (bwt.size() != bwtRev.size()) {
            throw std::runtime_error("bwt don't have the same size");
        }
        if (occ.size() != occRev.size()) {
            throw std::runtime_error("occ don't have the same size");
        }
        for (size_t sym{0}; sym < Sigma; ++sym) {
            if (occ.rank(occ.size(), sym) != occRev.rank(occ.size(), sym)) {
                throw std::runtime_error("wrong rank");
            }
        }
    }

    size_t memoryUsage() const {
        return occ.memoryUsage() + occRev.memoryUsage() + csa.memoryUsage();
    }

    size_t size() const {
        return occ.size();
    }

    size_t locate(size_t idx) const {
        auto opt = csa.value(idx);
        uint64_t steps{};
        while(!opt) {
            idx = occ.rank(idx, occ.symbol(idx));
            steps += 1;
            opt = csa.value(idx);
        }
        return opt.value() + steps;
    }
};

template <typename Index>
struct BiFMIndexCursor {
    static size_t constexpr Sigma = Index::Sigma;

    Index const* index;
    size_t lb;
    size_t lbRev;
    size_t len;
    BiFMIndexCursor()
        : index{nullptr}
    {}
    BiFMIndexCursor(Index const& index)
        : BiFMIndexCursor{index, 0, 0, index.size()}
    {}
    BiFMIndexCursor(Index const& index, size_t lb, size_t lbRev, size_t len)
        : index{&index}
        , lb{lb}
        , lbRev{lbRev}
        , len{len}
    {}
    bool empty() const {
        return len == 0;
    }
    size_t count() const {
        return len;
    }
    auto extendLeft() const -> std::array<BiFMIndexCursor, Sigma> {
//        prefetchLeft();
        auto const& occ = index->occ;
        occ.prefetch(lb+len);
        auto [rs1, prs1] = occ.all_ranks(lb);
        auto [rs2, prs2] = occ.all_ranks(lb+len);

        auto cursors = std::array<BiFMIndexCursor, Sigma>{};
        cursors[0] = BiFMIndexCursor{*index, rs1[0], lbRev, rs2[0] - rs1[0]};
        for (size_t i{1}; i < Sigma; ++i) {
            cursors[i] = BiFMIndexCursor{*index, rs1[i], lbRev + prs2[i-1] - prs1[i-1], rs2[i] - rs1[i]};
            //cursors[i].prefetchLeft();
        }
        cursors[1].prefetchLeft();
        return cursors;
    }

    auto extendRight() const -> std::array<BiFMIndexCursor, Sigma> {
//        prefetchRight();
        auto const& occ = index->occRev;
        occ.prefetch(lbRev+len);
        auto [rs1, prs1] = occ.all_ranks(lbRev);
        auto [rs2, prs2] = occ.all_ranks(lbRev+len);

        auto cursors = std::array<BiFMIndexCursor, Sigma>{};
        cursors[0] = BiFMIndexCursor{*index, lb, rs1[0], rs2[0] - rs1[0]};
        for (size_t i{1}; i < Sigma; ++i) {
            cursors[i] = BiFMIndexCursor{*index, lb + prs2[i-1] - prs1[i-1], rs1[i], rs2[i] - rs1[i]};
            //cursors[i].prefetchRight();
        }
        cursors[1].prefetchRight();
        return cursors;
    }
    void prefetchLeft() const {
        auto& occ = index->occ;
        occ.prefetch(lb);
        occ.prefetch(lb+len);
    }
    void prefetchRight() const {
        auto& occ = index->occRev;
        occ.prefetch(lbRev);
        occ.prefetch(lbRev+len);

    }

    auto extendLeft(uint8_t symb) const -> BiFMIndexCursor {
        assert(symb > 0);
//        prefetchLeft();
        auto& occ = index->occ;
        //occ.prefetch(lb+len);
        size_t newLb    = occ.rank(lb, symb);
        size_t newLbRev = lbRev + occ.prefix_rank(lb+len, symb-1) - occ.prefix_rank(lb, symb-1);
        size_t newLen   = occ.rank(lb+len, symb) - newLb;
        auto newCursor = BiFMIndexCursor{*index, newLb, newLbRev, newLen};
        newCursor.prefetchLeft();
        return newCursor;
    }
    auto extendRight(uint8_t symb) const -> BiFMIndexCursor {
        assert(symb > 0);
//        prefetchRight();
        auto& occ = index->occRev;
        //occ.prefetch(lbRev+len);
        size_t newLb    = lb + occ.prefix_rank(lbRev+len, symb-1) - occ.prefix_rank(lbRev, symb-1);
        size_t newLbRev = occ.rank(lbRev, symb);
        size_t newLen   = occ.rank(lbRev+len, symb) - newLbRev;
        auto newCursor = BiFMIndexCursor{*index, newLb, newLbRev, newLen};
        newCursor.prefetchRight();
        return newCursor;
    }

};
