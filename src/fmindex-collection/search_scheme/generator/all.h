// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "backtracking.h"
#include "bestKnown.h"
#include "h2.h"
#include "hato.h"
#include "kianfar.h"
#include "kucherov.h"
#include "lam.h"
#include "optimum.h"
#include "pex-bottom-up.h"
#include "pex-top-down.h"
#include "pigeon.h"
#include "suffixFilter.h"
#include "zeroOnesZero.h"

#include <functional>
#include <string>
#include <tuple>
#include <map>
#include <unordered_map>
#include <vector>

namespace fmc::search_scheme::generator {

struct GeneratorEntry {
    std::string name;
    std::string description;
    std::function<Scheme(int minError, int maxError, int sigma, int refSize)> generator;
};

inline auto all = []() {
    auto res = std::map<std::string, GeneratorEntry>{};
    auto add = [&](GeneratorEntry entry) {
        res.try_emplace(entry.name, entry);
    };
    add({ .name        = "backtracking",
          .description = "simple backtracking, not utilisying the bidirectional fm-index or search schemes",
          .generator   = []([[maybe_unused]] int minError, [[maybe_unused]] int maxError, [[maybe_unused]] int sigma, [[maybe_unused]] int dbSize) { return backtracking(1, minError, maxError); }
    });
    add({ .name        = "optimum",
          .description = "known optimim search schemes",
          .generator   = []([[maybe_unused]] int minError, [[maybe_unused]] int maxError, [[maybe_unused]] int sigma, [[maybe_unused]] int dbSize) { return optimum(minError, maxError); }
    });
    add({ .name        = "01*0",
          .description = "based on 01*0 seeds",
          .generator   = []([[maybe_unused]] int minError, [[maybe_unused]] int maxError, [[maybe_unused]] int sigma, [[maybe_unused]] int dbSize) { return zeroOnesZero_trivial(minError, maxError); }
    });
    add({ .name        = "01*0_opt",
          .description = "based on 01*0 seeds, but joining search schemes with same part order",
          .generator   = []([[maybe_unused]] int minError, [[maybe_unused]] int maxError, [[maybe_unused]] int sigma, [[maybe_unused]] int dbSize) { return zeroOnesZero_opt(minError, maxError); }
    });
    add({ .name        = "pigeon",
          .description = "based on the pigeon hole principle",
          .generator   = []([[maybe_unused]] int minError, [[maybe_unused]] int maxError, [[maybe_unused]] int sigma, [[maybe_unused]] int dbSize) { return pigeon_trivial(minError, maxError); }
    });
    add({ .name        = "pigeon_opt",
          .description = "based on the pigeon hole principle, removing duplicate paths",
          .generator   = []([[maybe_unused]] int minError, [[maybe_unused]] int maxError, [[maybe_unused]] int sigma, [[maybe_unused]] int dbSize) { return pigeon_opt(minError, maxError); }
    });
    add({ .name        = "suffix",
          .description = "based on suffix filter",
          .generator   = []([[maybe_unused]] int minError, [[maybe_unused]] int maxError, [[maybe_unused]] int sigma, [[maybe_unused]] int dbSize) { return suffixFilter(maxError+1, minError, maxError); }
    });
    add({ .name        = "lam",
          .description = "designed by kucherov, firstly designed by lam",
          .generator   = []([[maybe_unused]] int minError, [[maybe_unused]] int maxError, [[maybe_unused]] int sigma, [[maybe_unused]] int dbSize) { return lam(maxError); }
    });
    add({ .name        = "kianfar",
          .description = "designed by kianfar (?)",
          .generator   = []([[maybe_unused]] int minError, [[maybe_unused]] int maxError, [[maybe_unused]] int sigma, [[maybe_unused]] int dbSize) { return kianfar(maxError); }
    });
    add({ .name        = "kucherov-k1",
          .description = "designed by kucherov, divided into k+1 pieces (?)",
          .generator   = []([[maybe_unused]] int minError, [[maybe_unused]] int maxError, [[maybe_unused]] int sigma, [[maybe_unused]] int dbSize) { return kucherov(maxError+1, maxError); }
    });
    add({ .name        = "kucherov-k2",
          .description = "designed by kucherov, divided into k+2 pieces (?)",
          .generator   = []([[maybe_unused]] int minError, [[maybe_unused]] int maxError, [[maybe_unused]] int sigma, [[maybe_unused]] int dbSize) { return kucherov(maxError+2, maxError); }
    });
    add({ .name        = "h2-k1",
          .description = "designed by gottlieb, divided into k+1 pieces",
          .generator   = []([[maybe_unused]] int minError, [[maybe_unused]] int maxError, [[maybe_unused]] int sigma, [[maybe_unused]] int dbSize) { return h2(maxError+1, minError, maxError); }
    });
    add({ .name        = "h2-k2",
          .description = "designed by gottlieb, divided into k+2 pieces",
          .generator   = []([[maybe_unused]] int minError, [[maybe_unused]] int maxError, [[maybe_unused]] int sigma, [[maybe_unused]] int dbSize) { return h2(maxError+2, minError, maxError); }
    });
    add({ .name        = "h2-k3",
          .description = "designed by gottlieb, divided into k+3 pieces",
          .generator   = []([[maybe_unused]] int minError, [[maybe_unused]] int maxError, [[maybe_unused]] int sigma, [[maybe_unused]] int dbSize) { return h2(maxError+3, minError, maxError); }
    });
    add({ .name        = "hato",
          .description = "designed by hato",
          .generator   = []([[maybe_unused]] int minError, [[maybe_unused]] int maxError, [[maybe_unused]] int sigma, [[maybe_unused]] int dbSize) { return hato(maxError); }
    });
    add({ .name = "pex-bu",
          .description = "creates a search scheme, that mimics the pex bottom up strategy",
          .generator = []([[maybe_unused]] int minError, [[maybe_unused]] int maxError, [[maybe_unused]] int sigma, [[maybe_unused]] int dbSize) { return pex_bu(minError, maxError, false); }
    });
    add({ .name = "pex-td",
          .description = "creates a search scheme, that mimics the pex top down strategy",
          .generator = []([[maybe_unused]] int minError, [[maybe_unused]] int maxError, [[maybe_unused]] int sigma, [[maybe_unused]] int dbSize) { return pex_td(minError, maxError, false); }
    });
    add({ .name = "pex-bu-l",
          .description = "creates a search scheme, that mimics the pex bottom up strategy",
          .generator = []([[maybe_unused]] int minError, [[maybe_unused]] int maxError, [[maybe_unused]] int sigma, [[maybe_unused]] int dbSize) { return pex_bu(minError, maxError, true); }
    });
    add({ .name = "pex-td-l",
          .description = "creates a search scheme, that mimics the pex top down strategy",
          .generator = []([[maybe_unused]] int minError, [[maybe_unused]] int maxError, [[maybe_unused]] int sigma, [[maybe_unused]] int dbSize) { return pex_td(minError, maxError, true); }
    });

    return res;
}();

}
