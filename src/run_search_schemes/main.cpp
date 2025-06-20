// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <fmindex-collection/search_scheme/all.h>
#include <fmt/format.h>

namespace ss = fmc::search_scheme;

int main(int argc, char** argv) {
    if (argc != 4) {
        fmt::print("generators: ");
        for (auto const& [key, val] : ss::generator::all) {
            fmt::print("{}, ", key);
        }
        fmt::print("\n");
        fmt::print("call:\n");
        fmt::print("{} <len> <K> <gen>\n", argv[0]);
        return 0;
    }

    auto len = std::stod(argv[1]);
    auto K   = std::stod(argv[2]);
    auto gen = argv[3];
    {
        auto oss = ss::generator::all.at(gen).generator(0, K, 0, 0);
        if (oss.size() == 0) return 0;

//        for (auto s : oss) {
//            fmt::print("{}\n", fmt::join(s.pi, " "));
//            fmt::print("{}\n", fmt::join(s.l, " "));
//            fmt::print("{}\n\n", fmt::join(s.u, " "));
//        }
//        fmt::print("\n");

        auto ss = ss::expand(oss, len);

        auto nc = [&](auto ss) {
            return ss::weightedNodeCount</*Edit=*/false>(ss, 4, 3'000'000'000);
        };
        auto nce = [&](auto ss) {
            return ss::weightedNodeCount</*Edit=*/true>(ss, 4, 3'000'000'000);
        };
        auto dss = ss::expandByWNC</*Edit=*/true>(oss, len, 4, 3'000'000'000);


//        fmt::print("ess:\n");
//        for (auto s : ss) {
//            fmt::print("{}\n", fmt::join(s.pi, " "));
//            fmt::print("{}\n", fmt::join(s.l, " "));
//            fmt::print("{}\n\n", fmt::join(s.u, " "));
//        }
//        fmt::print("\ndss:\n");
//        for (auto s : dss) {
//            fmt::print("{}\n", fmt::join(s.pi, " "));
//            fmt::print("{}\n", fmt::join(s.l, " "));
//            fmt::print("{}\n\n", fmt::join(s.u, " "));
//        }


        auto valid = [&](auto ss) {
            return isValid(ss)/* and isComplete(ss, 0, K)*/;
        };

        fmt::print ("{:14}: {} ham/edit {:10.0f}/{:10.0f} → {:10.0f}/{:10.0f} ({}, {})\n", gen, K, nc(ss), nce(ss), nc(dss), nce(dss), valid(ss), valid(dss));
    }
}
