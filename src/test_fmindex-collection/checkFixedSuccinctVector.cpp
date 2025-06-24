// SPDX-FileCopyrightText: 2025, Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#include <catch2/catch_all.hpp>
#include <fstream>

#include "FixedSuccinctVector.h"
TEST_CASE("test FixedSuccinctVector", "[FixedSuccinctVector]") {
    SECTION("case1") {
        using Vec = fmc::FixedSuccinctVector<10>;

        auto vec = Vec{};
        vec.push_back(512);
        vec.push_back(17);
        vec.push_back(99);
        vec.push_back(0);
        vec.push_back(1023);
        vec.push_back(732);
        vec.push_back(101);

        auto expected = std::vector<size_t>{512, 17, 99, 0, 1023, 732, 101};
        auto check = [&]() {
            for (size_t i{0}; i < expected.size(); ++i) {
                INFO(i);
                CHECK(vec.at(i) == expected.at(i));
            }
        };

        check();

        auto change = [&](size_t i, size_t v) {
            vec.set(i, v);
            expected[i] = v;
            check();
        };
        change(0, 19);
        change(1, 20);
        change(2, 21);
        change(3, 22);
        change(4, 23);
        change(5, 24);
        change(6, 25);
    }
}
