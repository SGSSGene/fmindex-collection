#pragma once

#include <array>

namespace search_ng14 {

enum class Dir : uint8_t { Left, Right };
template <typename T>
struct Block {
    T       rank;
    size_t l;
    size_t u;
    Dir dir;
};

template <typename index_t, typename search_scheme_t, typename delegate_t>
struct Search {
    constexpr static size_t Sigma = index_t::Sigma;

    using cursor_t = BiFMIndexCursor<index_t>;
    using BlockIter = typename search_scheme_t::const_iterator;

    index_t const& index;
    search_scheme_t const& search;
    size_t qidx;
    delegate_t const& delegate;

    Search(index_t const& _index, search_scheme_t const& _search, size_t _qidx, delegate_t const& _delegate)
        : index     {_index}
        , search    {_search}
        , qidx      {_qidx}
        , delegate  {_delegate}
    {
        auto cur       = cursor_t{index};
        auto blockIter = search.begin();

        search_next<'M', 'M'>(cur, 0, blockIter);
    }

    template <bool Right>
    static auto extend(cursor_t const& cur, uint8_t symb) noexcept {
        if constexpr (Right) {
            return cur.extendRight(symb);
        } else {
            return cur.extendLeft(symb);
        }
    }

/*    template <bool Right, typename CB>
    static auto extend_cb(cursor_t const& cur, CB const& cb) noexcept {
        if constexpr (Right) {
            return cur.extendRight_cb(cb);
        } else {
            return cur.extendLeft_cb(cb);
        }
    }*/

    template <char LInfo, char RInfo>
    void search_next(cursor_t const& cur, int e, BlockIter blockIter) const noexcept {
        if (cur.count() == 0) {
            return;
        }

        if (blockIter == end(search)) {
            //if constexpr ((LInfo == 'M' or LInfo == 'I') and (RInfo == 'M' or RInfo == 'I')) {
                delegate(qidx, cur);
            //}
            /*if ((blockIter-1)->dir == Dir::Right) {
                search_next_dir_final<LInfo, RInfo, true>(cur, e, blockIter);
            } else {
                search_next_dir_final<LInfo, RInfo, false>(cur, e, blockIter);
            }*/
            return;
        }
        if (blockIter->dir == Dir::Right) {
            search_next_dir<LInfo, RInfo, true>(cur, e, blockIter);
        } else {
            search_next_dir<LInfo, RInfo, false>(cur, e, blockIter);
        }
    }

    template <char LInfo, char RInfo, bool Right>
    void search_next_dir(cursor_t const& cur, int e, BlockIter blockIter) const noexcept {
        static constexpr char TInfo = Right ? RInfo : LInfo;

        constexpr bool Deletion     = TInfo == 'M' or TInfo == 'D';
        constexpr bool Insertion    = TInfo == 'M' or TInfo == 'I';

        constexpr char OnMatchL      = Right ? LInfo : 'M';
        constexpr char OnMatchR      = Right ? 'M'   : RInfo;
        constexpr char OnSubstituteL = Right ? LInfo : 'S';
        constexpr char OnSubstituteR = Right ? 'S'   : RInfo;
        constexpr char OnDeletionL   = Right ? LInfo : 'D';
        constexpr char OnDeletionR   = Right ? 'D'   : RInfo;
        constexpr char OnInsertionL  = Right ? LInfo : 'I';
        constexpr char OnInsertionR  = Right ? 'I'   : RInfo;

        bool matchAllowed    = blockIter->l <= e and e <= blockIter->u;
        bool mismatchAllowed = blockIter->l <= e+1 and e+1 <= blockIter->u;

        auto symb = blockIter->rank;


        if (mismatchAllowed) {
            /*auto cursors = std::array<cursor_t, Sigma>{};
            for (size_t i{1}; i < Sigma; ++i) {
                if constexpr (Right) {
                    cursors[i] = cur.extendRight(i);
                } else {
                    cursors[i] = cur.extendLeft(i);
                }
            }*/
            auto cursors = [&]() {
                if constexpr (Right) {
                    return cur.extendRight();
                } else {
                    return cur.extendLeft();
                }
            }();

            /*extend_cb<Right>(cur, [&](auto c, auto newCur) {
                if (c > Sigma) return;
                cursors[c] = newCur;
            });*/
/*            auto cursors = [&]() {
                if constexpr (Right) {
                    return cur.extendRight();
                } else {
                    return cur.extendLeft();
               }
            }();*/


            if (matchAllowed) {
                auto newCur = cursors[symb];
                search_next<OnMatchL, OnMatchR>(newCur, e, blockIter+1);
            }

            for (uint8_t i{1}; i < symb; ++i) {
                auto newCur = cursors[i];

                if constexpr (Deletion) {
                    search_next<OnDeletionL, OnDeletionR>(newCur, e+1, blockIter); // deletion occurred in query
                }
                search_next<OnSubstituteL, OnSubstituteR>(newCur, e+1, blockIter+1); // as substitution
            }

            for (uint8_t i(symb+1); i < Sigma; ++i) {
                auto newCur = cursors[i];

                if constexpr (Deletion) {
                    search_next<OnDeletionL, OnDeletionR>(newCur, e+1, blockIter); // deletion occurred in query
                }
                search_next<OnSubstituteL, OnSubstituteR>(newCur, e+1, blockIter+1); // as substitution
            }


            if constexpr (Insertion) {
                //search_next_ins<OnInsertionL, OnInsertionR, Right>(cur, e+1, blockIter+1, cursors); // insertion occurred in query
                search_next<OnInsertionL, OnInsertionR>(cur, e+1, blockIter+1); // insertion occurred in query

            }


        } else if (matchAllowed) {
            auto newCur = extend<Right>(cur, symb);
            search_next<OnMatchL, OnMatchR>(newCur, e, blockIter+1);
        }
    }

    template <char LInfo, char RInfo, bool Right>
    void search_next_dir_final(cursor_t const& cur, int e, BlockIter blockIter) const noexcept {
        static constexpr char TInfo = Right ? RInfo : LInfo;

        bool mismatchAllowed = (blockIter-1)->l <= e+1 and e+1 <= (blockIter-1)->u;
        if (mismatchAllowed) {
            auto cursors = std::array<cursor_t, Sigma>{};
            extend_cb<Right>(cur, [&](auto c, auto newCur) {
                cursors[c] = newCur;
            });
            for (uint8_t i{1}; i < Sigma; ++i) {
                auto newCur = cursors[i];
                search_next<LInfo, RInfo>(newCur, e+1, blockIter); // deletion occurred in query
            }
        }
    }


    template <char LInfo, char RInfo, bool lastDir>
    void search_next_ins(cursor_t const& cur, int e, BlockIter blockIter, std::array<cursor_t, Sigma> const& cursors) const noexcept {
        if (cur.count() == 0) {
            return;
        }

        if (blockIter == end(search)) {
            if constexpr ((LInfo == 'M' or LInfo == 'I') and (RInfo == 'M' or RInfo == 'I')) {
                delegate(qidx, cur);
            }
            /*if ((blockIter-1)->dir == Dir::Right) {
                search_next_dir_final<LInfo, RInfo, true>(cur, e, blockIter);
            } else {
                search_next_dir_final<LInfo, RInfo, false>(cur, e, blockIter);
            }*/
            return;
        }
        if (lastDir == (blockIter->dir == Dir::Right)) {
            if (blockIter->dir == Dir::Right) {
                search_next_dir_ins<LInfo, RInfo, true>(cur, e, blockIter, cursors);
            } else {
                search_next_dir_ins<LInfo, RInfo, false>(cur, e, blockIter, cursors);
            }
        } else {
            if (blockIter->dir == Dir::Right) {
                search_next_dir<LInfo, RInfo, true>(cur, e, blockIter);
            } else {
                search_next_dir<LInfo, RInfo, false>(cur, e, blockIter);
            }

        }
    }

    template <char LInfo, char RInfo, bool Right>
    void search_next_dir_ins(cursor_t const& cur, int e, BlockIter blockIter, std::array<cursor_t, Sigma> const& cursors) const noexcept {
        static constexpr char TInfo = Right ? RInfo : LInfo;

        constexpr bool Deletion     = TInfo == 'M' or TInfo == 'D';
        constexpr bool Insertion    = TInfo == 'M' or TInfo == 'I';

        constexpr char OnMatchL      = Right ? LInfo : 'M';
        constexpr char OnMatchR      = Right ? 'M'   : RInfo;
        constexpr char OnSubstituteL = Right ? LInfo : 'S';
        constexpr char OnSubstituteR = Right ? 'S'   : RInfo;
        constexpr char OnDeletionL   = Right ? LInfo : 'D';
        constexpr char OnDeletionR   = Right ? 'D'   : RInfo;
        constexpr char OnInsertionL  = Right ? LInfo : 'I';
        constexpr char OnInsertionR  = Right ? 'I'   : RInfo;

        auto rank = blockIter->rank;

        bool matchAllowed    = blockIter->l <= e and e <= blockIter->u;
        bool mismatchAllowed = blockIter->l <= e+1 and e+1 <= blockIter->u;

        if (matchAllowed) {
            auto newCur = cursors[rank];
            search_next<OnMatchL, OnMatchR>(newCur, e, blockIter+1);
        }

        if (mismatchAllowed) {
            for (uint8_t i{1}; i < Sigma; ++i) {
                auto newCur = cursors[i];

                if constexpr (Deletion) {
                    search_next<OnDeletionL, OnDeletionR>(newCur, e+1, blockIter); // deletion occurred in query
                }
                search_next<OnSubstituteL, OnSubstituteR>(newCur, e+1, blockIter+1); // as substitution
            }

            for (uint8_t i(rank+1); i < Sigma; ++i) {
                auto newCur = cursors[i];

                if constexpr (Deletion) {
                    search_next<OnDeletionL, OnDeletionR>(newCur, e+1, blockIter); // deletion occurred in query
                }
                search_next<OnSubstituteL, OnSubstituteR>(newCur, e+1, blockIter+1); // as substitution
            }


            if constexpr (Insertion) {
                search_next_ins<OnInsertionL, OnInsertionR, Right>(cur, e+1, blockIter+1, cursors); // insertion occurred in query
            }
        }
    }

};


template <typename index_t, typename queries_t, typename search_schemes_t, typename delegate_t>
void search(index_t const & index, queries_t && queries, search_schemes_t const & search_scheme, delegate_t && delegate)
{
    if (search_scheme.empty()) return;
    auto internal_delegate = [&delegate] (size_t qidx, auto const & it) {
        delegate(qidx, it);
    };

    std::vector<std::vector<Block<size_t>>> search_scheme2;
    for (auto s : search_scheme) {
        std::vector<Block<size_t>> search2;
        for (size_t i{0}; i < s.pi.size(); ++i) {
            auto dir = [&]() {
                if (i == 0) {
                    return s.pi[i] < s.pi[i+1]?Dir::Right:Dir::Left;
                } else {
                    return s.pi[i-1] < s.pi[i]?Dir::Right:Dir::Left;
                }
            }();
            search2.emplace_back(Block<size_t>{{}, s.l[i], s.u[i], dir});
        }
        search_scheme2.emplace_back(move(search2));
    }

    for (size_t i{0}; i < queries.size(); ++i) {
        auto const& query = queries[i];
        for (size_t j{0}; j < search_scheme.size(); ++j) {
            auto& search = search_scheme2[j];
            for (size_t k {0}; k < search.size(); ++k) {
                search[k].rank = query[search_scheme[j].pi[k]];
            }
            Search<std::decay_t<decltype(index)>, std::decay_t<decltype(search)>, std::decay_t<decltype(internal_delegate)>>{index, search, i, internal_delegate};
        }
    }

}

//!\}

} // namespace seqan3
