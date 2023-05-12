// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include "backtracking.h"
#include "bestKnown.h"
#include "h2.h"
#include "kianfar.h"
#include "kucherov.h"
#include "optimum.h"
#include "pigeon.h"
#include "suffixFilter.h"
#include "zeroOnesZero.h"

#include <functional>
#include <string>
#include <tuple>
#include <vector>
#include <map>

namespace search_schemes::generator {

inline auto all = std::map<std::string, std::function<Scheme(int, int, int, int)>>{
    { "backtracking",   []([[maybe_unused]] int minError, [[maybe_unused]] int maxError, [[maybe_unused]] int sigma, [[maybe_unused]] int dbSize) { return backtracking(1, minError, maxError); }},
    { "optimum",        []([[maybe_unused]] int minError, [[maybe_unused]] int maxError, [[maybe_unused]] int sigma, [[maybe_unused]] int dbSize) { return optimum(minError, maxError); }},
    { "01*0",           []([[maybe_unused]] int minError, [[maybe_unused]] int maxError, [[maybe_unused]] int sigma, [[maybe_unused]] int dbSize) { return zeroOnesZero_trivial(minError, maxError); }},
    { "01*0_opt",       []([[maybe_unused]] int minError, [[maybe_unused]] int maxError, [[maybe_unused]] int sigma, [[maybe_unused]] int dbSize) { return zeroOnesZero_opt(minError, maxError); }},
    { "pigeon",         []([[maybe_unused]] int minError, [[maybe_unused]] int maxError, [[maybe_unused]] int sigma, [[maybe_unused]] int dbSize) { return pigeon_trivial(minError, maxError); }},
    { "pigeon_opt",     []([[maybe_unused]] int minError, [[maybe_unused]] int maxError, [[maybe_unused]] int sigma, [[maybe_unused]] int dbSize) { return pigeon_opt(minError, maxError); }},
    { "suffix",         []([[maybe_unused]] int minError, [[maybe_unused]] int maxError, [[maybe_unused]] int sigma, [[maybe_unused]] int dbSize) { return suffixFilter(maxError+1, minError, maxError); }},
    { "kianfar",        []([[maybe_unused]] int minError, [[maybe_unused]] int maxError, [[maybe_unused]] int sigma, [[maybe_unused]] int dbSize) { return kianfar(maxError); }},
    { "kucherov-k1",    []([[maybe_unused]] int minError, [[maybe_unused]] int maxError, [[maybe_unused]] int sigma, [[maybe_unused]] int dbSize) { return kucherov(maxError+1, maxError); }},
    { "kucherov-k2",    []([[maybe_unused]] int minError, [[maybe_unused]] int maxError, [[maybe_unused]] int sigma, [[maybe_unused]] int dbSize) { return kucherov(maxError+2, maxError); }},
    { "h2-k1",          []([[maybe_unused]] int minError, [[maybe_unused]] int maxError, [[maybe_unused]] int sigma, [[maybe_unused]] int dbSize) { return h2(maxError+1, minError, maxError); }},
    { "h2-k2",          []([[maybe_unused]] int minError, [[maybe_unused]] int maxError, [[maybe_unused]] int sigma, [[maybe_unused]] int dbSize) { return h2(maxError+2, minError, maxError); }},
    { "h2-k3",          []([[maybe_unused]] int minError, [[maybe_unused]] int maxError, [[maybe_unused]] int sigma, [[maybe_unused]] int dbSize) { return h2(maxError+3, minError, maxError); }},
};

}
