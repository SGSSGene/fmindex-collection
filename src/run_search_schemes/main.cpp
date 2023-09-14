// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#include <fmt/format.h>
#include <search_schemes/search_schemes.h>

int main(int argc, char** argv) {
    if (argc != 4) {
        fmt::print("generators: ");
        for (auto const& [key, val] : search_schemes::generator::all) {
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
        auto oss = search_schemes::generator::all.at(gen).generator(0, K, 0, 0);
        if (oss.size() == 0) return 0;

//        for (auto s : oss) {
//            fmt::print("{}\n", fmt::join(s.pi, " "));
//            fmt::print("{}\n", fmt::join(s.l, " "));
//            fmt::print("{}\n\n", fmt::join(s.u, " "));
//        }
//        fmt::print("\n");

        auto ss = search_schemes::expand(oss, len);

        auto nc = [&](auto ss) {
            return search_schemes::weightedNodeCount</*Edit=*/false>(ss, 4, 3'000'000'000);
        };
        auto nce = [&](auto ss) {
            return search_schemes::weightedNodeCount</*Edit=*/true>(ss, 4, 3'000'000'000);
        };
        auto dss = search_schemes::expandByWNC</*Edit=*/true>(oss, len, 4, 3'000'000'000);


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
