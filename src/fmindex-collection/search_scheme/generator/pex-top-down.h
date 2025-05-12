// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "pex-bottom-up.h"
#include "../isComplete.h"

namespace fmindex_collection::search_scheme::generator {

namespace pex {
inline auto constructTD(size_t maxError) -> std::unique_ptr<pex::Node> {
    auto root = std::make_unique<Node>();
    root->maxError = maxError;
    root->range = {0, maxError};

    auto openNodes = std::vector<Node*>{};
    openNodes.push_back(root.get());
    while (!openNodes.empty()) {
        auto top = openNodes.back();
        openNodes.pop_back();

        auto [start, end] = top->range;
        if (start == end) { // it is a leaf, return
            assert(top->maxError == 0);
            top->data = start;
            continue;
        }
        assert(top->maxError > 0);
        auto mid = (start+end)/2;

        auto lhs = std::make_unique<Node>();
        lhs->parent = top;
        lhs->maxError = top->maxError / 2;
        lhs->range = {start, mid};
        auto rhs = std::make_unique<Node>();
        rhs->parent = top;
        rhs->maxError = top->maxError - lhs->maxError;
        if (rhs->maxError > 0) rhs->maxError -= 1;
        rhs->range = {mid+1, end};

        openNodes.push_back(lhs.get());
        openNodes.push_back(rhs.get());

        auto children = Node::Children{};
        children.emplace_back(std::move(lhs));
        children.emplace_back(std::move(rhs));
        top->data = std::move(children);
    }
    return root;
};

}

inline auto pex_td(size_t minK, size_t K, bool increaseL) -> Scheme {
    [[maybe_unused]] auto N = K+1;

    assert(N>0);
    assert(minK <= K);
    assert(N>K);

    auto tree = pex::constructTD(K);
    tree->maxError = K; // max Error of the tree might be larger than required
    auto res = construct_scheme(*tree);

    // set minimum value
    for (auto& s : res) {
        s.l.back() = std::max(s.l.back(), minK);
    }

    // do some easy adjustemnts, by randomly increasing lower bound
    if (increaseL) {
        for (size_t i{0}; i < res.size(); ++i) {
            for (size_t j{res[i].l.size()}; j > 0; --j) {
                while (true) {
                    res[i].l[j-1] += 1;
                    if (!isComplete(res, minK, K)) {
                        res[i].l[j-1] -= 1;
                        break;
                    }
                }
            }
        }
    }

    return res;
}

}
