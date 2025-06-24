// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../Scheme.h"
#include "../isComplete.h"

#include <cassert>
#include <memory>
#include <tuple>
#include <variant>
#include <vector>

namespace fmc::search_scheme::generator {

namespace pex {
struct Node {
    Node* parent{nullptr};
    size_t maxError;
    std::tuple<size_t, size_t> range;

    using PartId   = size_t;
    using Children = std::vector<std::unique_ptr<Node>>;
    // `data` carries either a part number or a list of children
    std::variant<Children, PartId> data;

    void addChild(std::unique_ptr<Node> node) {
        node->parent = this;
        auto& children = std::get<Children>(data);
        if (children.empty()) {
            range = node->range;
        } else {
            auto  [_start, _end] = node->range;
            auto& [start, end]   = range;
            start = std::min(_start, start);
            end   = std::max(_end, end);
        }
        children.emplace_back(std::move(node));

        maxError = children.size()-1;
        for (auto& c : children) {
            maxError += c->maxError;
        }
    }

    void getAllLeafes(std::vector<Node const*>& retValue) const {
        if (std::get_if<Node::PartId>(&data)) {
            retValue.emplace_back(this);
            return;
        }

        for (auto& child : std::get<Node::Children>(data)) {
            child->getAllLeafes(retValue);
        }
    }

    auto getAllLeafes() const -> std::vector<Node const*> {
        auto retValue = std::vector<Node const*>{};
        getAllLeafes(retValue);
        return retValue;
    }
};

inline auto constructBU(size_t maxError) -> std::unique_ptr<Node> {
    size_t leafCt = maxError +1;
    auto nodes = std::vector<std::unique_ptr<Node>>{};
    for (size_t i{0}; i < leafCt; ++i) {
        nodes.emplace_back(std::make_unique<Node>());
        nodes.back()->maxError = 0;
        nodes.back()->data = i;
        nodes.back()->range = {i, i};
    }

    // merge nodes of level together, form binary with exception of last
    while (nodes.size() > 1) {
        auto newLevel = std::vector<std::unique_ptr<Node>>{};
        while (nodes.size() > 3) {
            auto node = std::make_unique<Node>();
            node->addChild(std::move(nodes[0]));
            node->addChild(std::move(nodes[1]));
            nodes.erase(nodes.begin(), nodes.begin()+2);
            newLevel.emplace_back(std::move(node));
        }
        // add the rest of nodes (might merge 3)
        {
            auto node = std::make_unique<Node>();
            for (auto& c : nodes) {
                node->addChild(std::move(c));
            }
            nodes.clear();
            newLevel.emplace_back(std::move(node));
        }
        std::swap(nodes, newLevel);
    }
    assert(nodes.size() == 1);
    auto node = std::unique_ptr<Node>{std::move(nodes[0])};
    return node;
};

}

inline auto construct_scheme(pex::Node const& tree) -> Scheme {
    auto leafs = tree.getAllLeafes();

    auto res = Scheme{};

    for (size_t i{0}; i < leafs.size(); ++i) {
        // generate search, starting with part i
        auto s = Search{};

        size_t minP{i+1}, maxP{i};
        auto ptr = leafs[i];
        while (ptr != nullptr) {
            auto [start, end] = ptr->range;

            // add leafs to the left
            if (start < minP) {
                for (size_t j{minP}; j > start; --j) {
                    s.pi.push_back(j-1);
                    s.l.push_back(0);
                    s.u.push_back(ptr->maxError);
                }
                minP = start;
            }
            // add leafs to the right
            if (end > maxP) {
                for (size_t j{maxP+1}; j <= end; ++j) {
                    s.pi.push_back(j);
                    s.l.push_back(0);
                    s.u.push_back(ptr->maxError);
                }
                maxP = end;
            }
            ptr = ptr->parent;
        }
        res.push_back(s);
    }
    return res;
}


inline auto pex_bu(size_t minK, size_t K, bool increaseL) -> Scheme {
    [[maybe_unused]] auto N = K+1;

    assert(N>0);
    assert(minK <= K);
    assert(N>K);

    auto tree = pex::constructBU(K);
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
