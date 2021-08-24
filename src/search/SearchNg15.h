#pragma once

#include <array>

namespace search_ng15 {

enum class Dir : uint8_t { Left, Right };

struct Block {
    size_t pi;
    size_t l;
    size_t u;
};

template <typename index_t, typename search_scheme_t, typename delegate_t>
struct Search {
    constexpr static size_t Sigma = index_t::Sigma;

    using cursor_t = BiFMIndexCursor<index_t>;
    using BlockIter = typename search_scheme_t::const_iterator;

    index_t const& index;
    search_scheme_t const& search;
    size_t qidx;
    std::vector<std::uint8_t> const& query;
    delegate_t const& delegate;
    size_t changePos;

    Search(index_t const& _index, search_scheme_t const& _search, size_t _qidx, std::vector<std::uint8_t> const& _query, delegate_t const& _delegate, size_t _changePos)
        : index     {_index}
        , search    {_search}
        , query     {_query}
        , qidx      {_qidx}
        , delegate  {_delegate}
        , changePos {_changePos}
    {
        auto cur       = cursor_t{index};
        auto blockIter = search.begin();

        if (changePos == 0) {
            search_next_left<'M'>(cur, 0, 0);
        } else {
            search_next_right<'M'>(cur, 0, 0);
        }
    }

    template <Dir Direction>
    static auto extend(cursor_t const& cur, uint8_t symb) noexcept {
        if constexpr (Direction == Dir::Right) {
            return cur.extendRight(symb);
        } else {
            return cur.extendLeft(symb);
        }
    }
    template <Dir Direction>
    static auto extend(cursor_t const& cur) noexcept {
        if constexpr (Direction == Dir::Right) {
            return cur.extendRight();
        } else {
            return cur.extendLeft();
        }
        /*auto cursors = std::array<cursor_t, Sigma>{};
        for (size_t i{1}; i < Sigma; ++i) {
            if constexpr (Direction == Dir::Right) {
                cursors[i] = cur.extendRight(i);
            } else {
                cursors[i] = cur.extendLeft(i);
            }
        }
        return cursors;*/
    }




    template <char Info>
    void search_next_right(cursor_t const& cur, size_t e, size_t pos) const noexcept {
        if (cur.count() == 0) {
            return;
        }

        if (pos == changePos) {
            if constexpr ((Info == 'M' or Info == 'S')) {
                search_next_left<'M'>(cur, e, pos);
            }
            return;
        }
        search_next_dir<Info, Dir::Right>(cur, e, pos);
    }

    template <char Info>
    void search_next_left(cursor_t const& cur, size_t e, size_t pos) const noexcept {
        if (cur.count() == 0) {
            return;
        }

        if (pos == query.size()) {
            if constexpr ((Info == 'M' or Info == 'S')) {
                delegate(qidx, cur, e);
            }
            return;
        }
        search_next_dir<Info, Dir::Left>(cur, e, pos);
    }

    template <char Info, Dir Direction>
    void search_next(cursor_t const& cur, size_t e, size_t pos) const noexcept {
        if constexpr (Direction == Dir::Right) {
            search_next_right<Info>(cur, e, pos);
        } else {
            search_next_left<Info>(cur, e, pos);
        }
    }


    template <char LastOp, Dir Direction>
    void search_next_dir(cursor_t const& cur, size_t e, size_t pos) const noexcept {

        constexpr bool Deletion     = LastOp == 'M' or LastOp == 'D';
        constexpr bool Insertion    = LastOp == 'M' or LastOp == 'I';

        auto const& block = search[pos];
        bool matchAllowed    = block.l <= e and e <= block.u;
        bool mismatchAllowed = block.l <= e+1 and e+1 <= block.u;

        auto symb = query[block.pi];

        if (mismatchAllowed) {
            auto cursors = extend<Direction>(cur);

            if (matchAllowed) {
                auto newCur = cursors[symb];
                search_next<'M', Direction>(newCur, e, pos+1);
            }

            for (uint8_t i{1}; i < symb; ++i) {
                auto newCur = cursors[i];

                if constexpr (Deletion) {
                    search_next<'D', Direction>(newCur, e+1, pos); // deletion occurred in query
                }
                search_next<'S', Direction>(newCur, e+1, pos+1); // as substitution
            }

            for (uint8_t i(symb+1); i < Sigma; ++i) {
                auto newCur = cursors[i];

                if constexpr (Deletion) {
                    search_next<'D', Direction>(newCur, e+1, pos); // deletion occurred in query
                }
                search_next<'S', Direction>(newCur, e+1, pos+1); // as substitution
            }


            if constexpr (Insertion) {
                search_next<'I', Direction>(cur, e+1, pos+1); // insertion occurred in query
            }

        } else if (matchAllowed) {
            auto newCur = extend<Direction>(cur, symb);
            search_next<'M', Direction>(newCur, e, pos+1);
        }
    }


};

template <typename index_t, typename queries_t, typename search_schemes_t, typename delegate_t>
void search(index_t const & index, queries_t && queries, search_schemes_t const & search_scheme, delegate_t && delegate)
{
    if (search_scheme.empty()) return;

    std::vector<size_t> dirChange;
    std::vector<std::vector<Block>> search_scheme2;
    for (auto s : search_scheme) {
        std::vector<Block> search2;
        for (size_t i{0}; i < s.pi.size(); ++i) {
            search2.emplace_back(Block{(size_t)s.pi[i], (size_t)s.l[i], (size_t)s.u[i]});
        }

        bool startDir = s.pi[0] < s.pi[1];
        size_t change = s.pi.size();
        for (size_t i{1}; i < s.pi.size(); ++i) {
            auto dir = s.pi[i-1] < s.pi[i];
            if (dir != startDir) {
                change = i;
                break;
            }
        }
        if (s.pi[0] > s.pi[1]) {
            change = 0;
        }
        dirChange.push_back(change);

        search_scheme2.emplace_back(move(search2));
    }

    for (size_t i{0}; i < queries.size(); ++i) {
        auto const& query = queries[i];
        for (size_t j{0}; j < search_scheme.size(); ++j) {
            auto& search = search_scheme2[j];
            Search{index, search, i, query, [&](size_t qidx, auto const& it, size_t e) {
                delegate(qidx, it, e);
            }, dirChange[j]};
        }
    }

}

}
