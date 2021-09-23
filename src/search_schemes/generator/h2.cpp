#include "h2.h"

#include <cassert>
#include <cstddef>

namespace oss::generator {
namespace {


auto pi(int row, int n, int N, int K, int Mod) {
    row = K - row;

    int shiftRight = Mod * row;
    n = n + shiftRight;

    if (n < N-row) {
        return n + row;
    }
    return N+shiftRight - n-1;
};

auto generatePieces(int N, int K, int Mod=0) {
    auto pieces = std::vector<std::vector<int>>(K+1, std::vector<int>(N, 0));
    for (int row{0}; row < K+1; ++row) {
        for (int i{0}; i < N; ++i) {
            pieces[row][i] = pi(row, i, N, K, Mod);
        }
    }
    return pieces;
}

auto generateDiffMatrix(int N, int K) {
    auto diffs = std::vector<std::vector<int>>(K+1, std::vector<int>(N, 0));
    for (int i{K}; i<N; ++i) {
        for (int row{0}; row < K+1; ++row) {
            diffs[row][i] = K-row;
        }
    }
    for (int i{0}; i<K; ++i) {
        for (int row{0}; row < K; ++row) {
            diffs[row][i] = (row-i + K) % K;
        }
        diffs[K][i] = K;

    }
    return diffs;
}

auto generateOptimizedDiffMatrix(int N, int K) {
    auto mat = generateDiffMatrix(N, K);

    auto isValid = [&](int row, int n, int v) {
        if (row == n) {
            return false;
        }
        if (row > n) {
            for (int i{0}; i<n; ++i) {
                if (mat[row][i] < v) {
                    return false;
                }
            }
        } else {
            for (int i{row+1}; i<n; ++i) {
                if (mat[row][i] > v) {
                    return false;
                }
            }
        }
        return true;
    };

    for (int i{0}; i<N; ++i) {
        for (int j{0}; j<K+1; ++j) {
            if (i == j or mat[j][i] == 0) continue;
            if (not isValid(j, i, mat[j][i])) {
                // check for switch
                int index=-1;
                for (int k{j+1}; k<K+1; ++k) {
                    if (isValid(j, i, mat[k][i]) and isValid(k, i, mat[j][i])) {
                        index = k;
                        break;
                    }
                }
                assert(index != -1);
                std::swap(mat[index][i], mat[j][i]);
            }
        }
    }

    return mat;
}


auto generateLowerBound(int N, int K) {
    auto bound = std::vector<std::vector<int>>(K+1, std::vector<int>(N, 0));
    for (int i{0}; i <= K; ++i) {
        for (int j{0}; j<K-i+1; ++j) {
            bound[i][N-j-1] = i;
        }
    }
    return bound;
}

auto generateUpperBound(std::vector<std::vector<int>> const& pieces, std::vector<std::vector<int>> const& lower) {
    assert(pieces.size() >= 1);
    int const K = pieces.size() - 1;
    int const N = pieces[0].size();
    assert(N >= K);

    auto diffs = generateOptimizedDiffMatrix(N, K);

    auto bound = std::vector<std::vector<int>>(K+1, std::vector<int>(N, 0));

    for (int i{1}; i<N; ++i) {
        for (int row{K}; row >= 0; --row) {
            assert(pieces[row].size() == N);
            int j = pieces[row][i];
            bound[row][i] = std::max(bound[row][i-1], lower[row][i-1] + diffs[K-row][j]);
        }
    }
    return bound;
}

}

auto h2(int N, int minK, int K) -> Scheme {
    assert(N>0);
    assert(minK >= 0);
    assert(minK <= K);
    assert(K>=0);
    assert(N>K);

    auto pieces = generatePieces(N, K, 0);
    auto lower  = generateLowerBound(N, K);
    auto upper  = generateUpperBound(pieces, lower);
    auto ss = Scheme{};
    for (size_t i{0}; i < pieces.size(); ++i) {
        auto s = SearchTree {pieces[i], lower[i], upper[i]};
        for (auto& pi : s.pi) {
            pi += 1;
        }
        ss.push_back(s);
    }

    // set minimum value
    for (auto& s : ss) {
        s.l.back() = std::max(s.l.back(), minK);
    }

//    assert(isValid(ss));

    return ss;
}
}
