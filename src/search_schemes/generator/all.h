#pragma once

#include "backtracking.h"
#include "bestKnown.h"
#include "h2.h"
#include "optimum.h"
#include "pigeon.h"
#include "suffixFilter.h"
#include "zeroOnesZero.h"
#include "kianfar.h"
#include "kucherov.h"

#include <functional>
#include <string>
#include <tuple>
#include <vector>
#include <map>

namespace oss::generator {

inline auto all = std::map<std::string, std::function<Scheme(int, int, int, int)>>{
    { "backtracking",   [](int minError, int maxError, int sigma, int dbSize) { return oss::generator::backtracking(1, minError, maxError); }},
    { "optimum",        [](int minError, int maxError, int sigma, int dbSize) { return oss::generator::optimum(maxError+1, minError, maxError); }},
    { "01*0",           [](int minError, int maxError, int sigma, int dbSize) { return oss::generator::zeroOnesZero_trivial(minError, maxError); }},
    { "01*0_opt",       [](int minError, int maxError, int sigma, int dbSize) { return oss::generator::zeroOnesZero_opt(minError, maxError); }},
    { "pigeon",         [](int minError, int maxError, int sigma, int dbSize) { return oss::generator::pigeon_trivial(minError, maxError); }},
    { "pigeon_opt",     [](int minError, int maxError, int sigma, int dbSize) { return oss::generator::pigeon_opt(minError, maxError); }},
    { "suffix",         [](int minError, int maxError, int sigma, int dbSize) { return oss::generator::suffixFilter(maxError+1, minError, maxError); }},
    { "kianfar",        [](int minError, int maxError, int sigma, int dbSize) { return oss::generator::kianfar(minError, maxError); }},
    { "kucherov-k1",    [](int minError, int maxError, int sigma, int dbSize) { return oss::generator::kucherov(maxError+1, minError, maxError); }},
    { "kucherov-k2",    [](int minError, int maxError, int sigma, int dbSize) { return oss::generator::kucherov(maxError+2, minError, maxError); }},
    { "h2-k1",          [](int minError, int maxError, int sigma, int dbSize) { return oss::generator::h2(maxError+1, minError, maxError); }},
    { "h2-k2",          [](int minError, int maxError, int sigma, int dbSize) { return oss::generator::h2(maxError+2, minError, maxError); }},
    { "h2-k3",          [](int minError, int maxError, int sigma, int dbSize) { return oss::generator::h2(maxError+3, minError, maxError); }},
};

}
