// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "FMIndex.h"

namespace fmc {

/** Special FMIndex, that
 * only carries entries with the same length.
 *
 * It always searches for the first character and then the next in a fixed order and no overlappings
 *
 */
template <size_t TSigma, template <size_t> typename String = string::FlattenedBitvectors_512_64k>
    requires String_c<String<TSigma>>
struct LinearFMIndex {
    static size_t constexpr Sigma = TSigma;

    struct Column {
        String<Sigma>               bwt;
        std::array<size_t, Sigma+1> C{0};
    };
    std::vector<Column> columns;

    size_t size_;
    size_t depth_;

    std::vector<size_t> ordered{};


    LinearFMIndex(Sequences auto const& _inputs) {
        auto lastOrder = std::vector<size_t>{};
        for (size_t i{0}; i < _inputs.size(); ++i) {
            lastOrder.push_back(i);
        }

        // counts characters in a column
        auto countColumn = [&](size_t col, size_t begin, size_t end) {
            auto count = std::array<size_t, 257>{};
            for (auto i{begin}; i < end; ++i) {
                auto row = lastOrder[i];
                auto c = _inputs[row][col];
                count[c] += 1;
            }
            return count;
        };
        auto accCount = [](std::array<size_t, 257>& count) {
            for (size_t i{1}; i < count.size(); ++i) {
                count[i] = count[i-1] + count[i];
            }
        };

        auto target = std::vector<std::vector<uint8_t>>{};
        target.resize(_inputs.size(), std::vector<uint8_t>(_inputs[0].size()));


        auto sort = std::function<void(size_t, size_t, size_t)>{};

        auto pos = std::vector<size_t>{};
        pos.resize(lastOrder.size());

        columns.resize(_inputs[0].size());

        auto temp_input = std::vector<uint8_t>{};
        temp_input.resize(_inputs.size());

        for (size_t j{0}; j < _inputs[0].size(); ++j) {
            size_t col = _inputs[0].size() - j - 1;
            auto count = countColumn(col, 0, _inputs.size());
            accCount(count);

            for (size_t _i{0}; _i < _inputs.size(); ++_i) {
                auto i = _inputs.size() - _i - 1;
                auto row = lastOrder[i];
                auto c = _inputs[row][col];
                pos[--count[c]] = row;
            }
            std::swap(lastOrder, pos);
            //fmt::print("pass {}\n", j);
            //printPos(lastOrder);


            // fill bwts + C
            auto tcol = ((col == 0)?columns.size():col)-1;
            auto& bwt = columns[tcol].bwt;
            auto& C = columns[tcol].C;

            // fill temp string for bwt
            for (size_t i{0}; i < lastOrder.size(); ++i) {
                auto row = lastOrder[i];
                auto c = _inputs[row].back();
                if (col > 0) c = _inputs[row][col-1];
                temp_input[i] = c;
            }
            bwt = {temp_input};

            // fill C
            for (size_t i{0}; i <= Sigma; ++i) {
                C[i] = bwt.prefix_rank(bwt.size(), i);
            }
        }
        size_ = _inputs.size();
        ordered = std::move(lastOrder);
    }

    auto size() const {
        return size_;
    }

    auto depth() const {
        return columns.size();
    }

    auto locate(size_t idx) const -> size_t {
        return ordered[idx];
    }
};

}
