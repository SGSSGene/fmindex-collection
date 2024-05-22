// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../Scheme.h"

#include <cassert>
#include <cstddef>
#include <limits>

namespace search_schemes::generator {

namespace h2_detail {

inline auto pi(size_t row, size_t n, size_t N, size_t K, size_t Mod) {
    row = K - row;

    size_t shiftRight = Mod * row;
    n = n + shiftRight;

    if (n < N-row) {
        return n + row;
    }
    return N+shiftRight - n-1;
}

inline auto generatePieces(size_t N, size_t K, size_t Mod=0) {
    auto pieces = std::vector<std::vector<size_t>>(K+1, std::vector<size_t>(N, 0));
    for (size_t row{0}; row < K+1; ++row) {
        for (size_t i{0}; i < N; ++i) {
            pieces[row][i] = pi(row, i, N, K, Mod);
        }
    }
    return pieces;
}

inline auto generateDiffMatrix(size_t N, size_t K) {
    auto diffs = std::vector<std::vector<size_t>>(K+1, std::vector<size_t>(N, 0));
    for (size_t i{K}; i<N; ++i) {
        for (size_t row{0}; row < K+1; ++row) {
            diffs[row][i] = K-row;
        }
    }
    for (size_t i{0}; i<K; ++i) {
        for (size_t row{0}; row < K; ++row) {
            diffs[row][i] = (row-i + K) % K;
        }
        diffs[K][i] = K;

    }
    return diffs;
}

inline auto generateOptimizedDiffMatrix(size_t N, size_t K) {
    auto mat = generateDiffMatrix(N, K);

    auto isValid = [&](size_t row, size_t n, size_t v) {
        if (row == n) {
            return false;
        }
        if (row > n) {
            for (size_t i{0}; i<n; ++i) {
                if (mat[row][i] < v) {
                    return false;
                }
            }
        } else {
            for (size_t i{row+1}; i<n; ++i) {
                if (mat[row][i] > v) {
                    return false;
                }
            }
        }
        return true;
    };

    for (size_t i{0}; i<N; ++i) {
        for (size_t j{0}; j<K+1; ++j) {
            if (i == j or mat[j][i] == 0) continue;
            if (not isValid(j, i, mat[j][i])) {
                // check for switch
                size_t index = std::numeric_limits<size_t>::max();
                for (size_t k{j+1}; k<K+1; ++k) {
                    if (isValid(j, i, mat[k][i]) and isValid(k, i, mat[j][i])) {
                        index = k;
                        break;
                    }
                }
                assert(index != std::numeric_limits<size_t>::max());
                std::swap(mat[index][i], mat[j][i]);
            }
        }
    }

    return mat;
}

inline auto generateLowerBound(size_t N, size_t K) {
    auto bound = std::vector<std::vector<size_t>>(K+1, std::vector<size_t>(N, 0));
    for (size_t i{0}; i <= K; ++i) {
        for (size_t j{0}; j<K-i+1; ++j) {
            bound[i][N-j-1] = i;
        }
    }
    return bound;
}

inline auto generateUpperBound(std::vector<std::vector<size_t>> const& pieces, std::vector<std::vector<size_t>> const& lower, size_t N, size_t K) {
    assert(pieces.size() >= 1);
    assert(N >= K);

    auto diffs = generateOptimizedDiffMatrix(N, K);

    auto bound = std::vector<std::vector<size_t>>(K+1, std::vector<size_t>(N, 0));

    for (size_t i{1}; i<N; ++i) {
        for (size_t _row{K+1}; _row > 0; --_row) {
            size_t row = _row - 1;
            assert(pieces[row].size() == N);
            size_t j = pieces[row][i];
            bound[row][i] = std::max(bound[row][i-1], lower[row][i-1] + diffs[K-row][j]);
        }
    }
    return bound;
}

inline auto h2(size_t N, size_t minK, size_t K) -> Scheme {
    assert(N>0);
    assert(minK <= K);
    assert(N >= K);

    auto pieces = generatePieces(N, K, 0);
    auto lower  = generateLowerBound(N, K);
    auto upper  = generateUpperBound(pieces, lower, N, K);
    auto ss = Scheme{};
    for (size_t i{0}; i < pieces.size(); ++i) {
        ss.emplace_back(Search{pieces[i], lower[i], upper[i]});
    }

    // set minimum value
    for (auto& s : ss) {
        s.l.back() = std::max(s.l.back(), minK);
    }

    return ss;
}

}

inline auto h2(size_t N, size_t minK, size_t K) -> Scheme {
    return h2_detail::h2(N, minK, K);
}

}
