// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "StopWatch.h"
#include "random.h"

#include <fmindex-collection/utils.h>

#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>



inline auto construct_bwt_from_sa(std::vector<uint64_t> const& sa, std::string_view const& text) -> std::vector<uint8_t> {
    assert(sa.size() == text.size());
    std::vector<uint8_t> bwt;
    bwt.resize(text.size());
    for (size_t i{0}; i < sa.size(); ++i) {
        bwt[i] = text[(sa[i] + text.size()- 1) % text.size()];
    }
    return bwt;
}
inline auto readFile(std::filesystem::path const& file) -> std::vector<uint8_t> {
    auto ifs = std::ifstream{file, std::ios::binary};
    ifs.seekg(0, std::ios::end);
    std::size_t size = ifs.tellg();
    auto buffer = std::vector<uint8_t>(size);
    ifs.seekg(0, std::ios::beg);
    ifs.read(reinterpret_cast<char*>(buffer.data()), buffer.size());
    return buffer;
}

template <size_t Sigma>
auto generateBWT(size_t len) -> std::tuple<std::vector<uint8_t>, std::vector<uint8_t>> {
    StopWatch watch;
    std::string text;
    text.resize(len, '\0');
    for (size_t i{0}; i < text.size(); ++i) {
        text[i] = (xorshf96() % (Sigma-1))+1;
    }
    text.back() = '\0';

    auto time_generation = watch.reset();
    std::cout << "text size: " << text.size()/1024/1024 << " million chars, "<<  text.size()/1024/1024*std::log(Sigma)/std::log(2) << " million bits\n";
    std::cout << "text generation: "<< time_generation << "s\n";

    auto bwt = [&]() {
        auto sa = fmc::createSA64({reinterpret_cast<uint8_t const*>(text.data()), text.size()}, 1);
        auto time_saconstruction = watch.reset();
        std::cout << "sa - construction time: "<< time_saconstruction << "s\n";

        return construct_bwt_from_sa(sa, text);
    }();
    auto bwtRev = [&]() {
        std::ranges::reverse(text);
        auto sa = fmc::createSA64({reinterpret_cast<uint8_t const*>(text.data()), text.size()}, 1);
        auto time_saconstruction = watch.reset();
        std::cout << "sa - rev construction time: "<< time_saconstruction << "s\n";

        return construct_bwt_from_sa(sa, text);
    }();
    return {bwt, bwtRev};
}


struct Result {
    std::string name;
    double expectedMemory                  = std::numeric_limits<double>::quiet_NaN();
    double constructionTime                = std::numeric_limits<double>::quiet_NaN();
    double benchV1                         = std::numeric_limits<double>::quiet_NaN();
    size_t benchV1CheckSum                 = 0;
    double benchV2                         = std::numeric_limits<double>::quiet_NaN();
    size_t benchV2CheckSum                 = 0;
    double benchV3                         = std::numeric_limits<double>::quiet_NaN();
    std::array<size_t, 2> benchV3CheckSum  = {0, 0};
    double benchV4                         = std::numeric_limits<double>::quiet_NaN();
    std::array<size_t, 2> benchV4CheckSum  = {0, 0};
    double benchV5                         = std::numeric_limits<double>::quiet_NaN();
    size_t benchV5CheckSum                 = 0;
    double benchV6                         = std::numeric_limits<double>::quiet_NaN();
    std::array<size_t, 2> benchV6CheckSum  = {0, 0};
    double totalTime                       = std::numeric_limits<double>::quiet_NaN();
};
