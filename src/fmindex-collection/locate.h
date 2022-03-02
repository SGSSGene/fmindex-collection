#pragma once

#pragma once

namespace fmindex_collection {


template <typename index_t, typename cursor_t>
struct LocateLinear {
    struct iter {
        index_t const& index;
        size_t pos;

        auto operator*() const -> std::tuple<size_t, size_t> {
            return index.locate(pos);
        }
        auto operator++() -> iter& {
            pos += 1;
            return *this;
        }
        auto operator!=(iter const& other) {
            return pos != other.pos;
        }
    };

    index_t const& index;
    cursor_t cursor;

    friend auto begin(LocateLinear const& locate) {
        return iter{locate.index, locate.cursor.lb};
    }
    friend auto end(LocateLinear const& locate) {
        return iter{locate.index, locate.cursor.lb + locate.cursor.len};
    }
};

template <typename index_t, typename cursor_t>
struct LocateFMTree {
    std::vector<std::tuple<size_t, size_t>> positions;

    LocateFMTree(index_t const& index, cursor_t cursor, size_t maxDepth) {
        positions.reserve(cursor.count());
        auto stack = std::vector<std::tuple<cursor_t, size_t>>{};
        stack.emplace_back(cursor, 0);
        size_t samplingRate = index.csa.samplingRate;
        while(!stack.empty()) {
            auto [cursor, depth] = stack.back();
            stack.pop_back();
            //!TODO what should the termination criteria be?
            if (depth >= maxDepth or cursor.count() < index_t::Sigma*2) {
                for (size_t pos{cursor.lb}; pos < cursor.lb + cursor.len; ++pos) {
                    auto v = index.locate(pos, samplingRate - depth - 1);
                    if (v) {
                        auto [seqid, seqpos] = *v;
                        positions.emplace_back(seqid, seqpos+depth);
                    }
                }
            } else {
                for (size_t pos{cursor.lb}; pos < cursor.lb + cursor.len; ++pos) {
                    auto v = index.single_locate_step(pos);
                    if (v) {
                        auto [seqid, seqpos] = *v;
                        positions.emplace_back(seqid, seqpos+depth);
                    }
                }

                auto cursors = cursor.extendLeft();
                for (size_t sym{1}; sym < cursors.size(); ++sym) {
                    stack.emplace_back(cursors[sym], depth+1);
                }
            }
        }
        assert(positions.size() == cursor.count());
    }

    friend auto begin(LocateFMTree const& locate) {
        return begin(locate.positions);
    }
    friend auto end(LocateFMTree const& locate) {
        return end(locate.positions);
    }
};


}
