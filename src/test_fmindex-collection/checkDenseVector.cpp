// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: CC0-1.0

#include <catch2/catch_all.hpp>
#include <fmindex-collection/DenseVector.h>
#include <fstream>

TEST_CASE("checking dense vector", "[densevector]") {
    auto input = std::vector<uint64_t>{10, 11, 5, 0, 6, 1, 7, 2, 3, 8, 4, 9};
    auto vec   = fmc::DenseVector{input};

    for (size_t i{0}; i < vec.size(); ++i) {
        CHECK(vec[i] == input[i]);
    }
}

TEST_CASE("checking dense vector with common divisor with small scaling factor", "[densevector]") {
    auto input = std::vector<uint64_t>{10, 12, 16, 0, 6, 22, 18, 2, 20, 8, 4, 10};
    auto vec   = fmc::DenseVector{input};

    for (size_t i{0}; i < vec.size(); ++i) {
        CHECK(vec[i] == input[i]);
    }
}

TEST_CASE("checking dense vector with common divisor with larger factor", "[densevector]") {
    auto input = std::vector<uint64_t>{100, 125, 160, 50, 60, 225, 180, 25, 200, 85, 40, 105};
    auto vec   = fmc::DenseVector{input};

    for (size_t i{0}; i < vec.size(); ++i) {
        CHECK(vec[i] == input[i]);
    }
}

TEST_CASE("checking dense vector with common divisor and push_back", "[densevector]") {
    auto input = std::vector<uint64_t>{100, 125, 160, 50, 60, 225, 180, 25, 200, 85, 40, 105};
    auto vec   = fmc::DenseVector(/*.largestValue =*/ 225, /*.commonDivisor = */ 5);
    for (auto i : input) {
        vec.push_back(i);
    }

    for (size_t i{0}; i < vec.size(); ++i) {
        CHECK(vec[i] == input[i]);
    }
}

TEST_CASE("checking dense vector concat operation", "[densevector]") {
    auto inputL = std::vector<uint64_t>{10, 12, 16, 0, 6, 22, 18, 2, 20, 8, 4, 10};
    auto vecL   = fmc::DenseVector{inputL};
    auto inputR = std::vector<uint64_t>{100, 125, 160, 50, 60, 225, 180, 25, 200, 85, 40, 105};
    auto vecR   = fmc::DenseVector{inputR};

    auto input = inputL;
    for (auto i : inputR) {
        input.push_back(i);
    }

    auto vec = fmc::DenseVector::concat(vecL, vecR);

    for (size_t i{0}; i < vec.size(); ++i) {
        CHECK(vec[i] == input[i]);
    }
}

TEST_CASE("checking dense vector concat operation with gcd", "[densevector]") {
    auto inputL = std::vector<uint64_t>{100, 120, 160, 0, 60, 220, 180, 20, 200, 80, 40, 100};
    auto vecL   = fmc::DenseVector{inputL};
    auto inputR = std::vector<uint64_t>{300, 375, 480, 150, 180, 675, 540, 75, 600, 255, 120, 315};
    auto vecR   = fmc::DenseVector{inputR};

    auto input = inputL;
    for (auto i : inputR) {
        input.push_back(i);
    }

    auto vec = fmc::DenseVector::concat(vecL, vecR);

    CHECK(vec.commonDivisor == 5);
    for (size_t i{0}; i < vec.size(); ++i) {
        CHECK(vec[i] == input[i]);
    }
}
