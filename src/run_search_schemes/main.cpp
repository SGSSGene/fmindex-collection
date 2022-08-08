#include <fmt/format.h>
#include <search_schemes/generator/all.h>
#include <search_schemes/nodeCount.h>
#include <search_schemes/expectedNodeCount.h>
#include <search_schemes/expand.h>
#include <search_schemes/isComplete.h>

int main(int argc, char** argv) {
    if (argc != 3) {
        fmt::print("generators: ");
        for (auto const& [key, val] : search_schemes::generator::all) {
            fmt::print("{}, ", key);
        }
        fmt::print("\n");
        return 0;
    }

    auto K   = std::stod(argv[1]);
    auto gen = argv[2];
    {
        auto oss = search_schemes::generator::all.at(gen)(0, K, 0, 0);
//        auto oss = search_schemes::generator::pigeon_opt(0, K);
/*        for (auto s : oss) {
            fmt::print("{}\n", fmt::join(s.pi, " "));
            fmt::print("{}\n", fmt::join(s.l, " "));
            fmt::print("{}\n\n", fmt::join(s.u, " "));
        }*/

        auto ss = search_schemes::expand(oss, 100);
//        fmt::print("node count ham: {}\n",  search_schemes::nodeCount(ss, 4));
//        fmt::print("node count edit: {}\n", search_schemes::nodeCountEdit(ss, 4));

        auto nc = [&](auto ss) {
            return search_schemes::expectedNodeCount(ss, 4, 3'000'000'000);
        };
        auto nce = [&](auto ss) {
            return search_schemes::expectedNodeCountEdit(ss, 4, 3'000'000'000);
        };
        auto ess = search_schemes::expandDynamic(oss, 100, 4, 3'000'000'000);

        auto valid = [&](auto ss) {
            return isValid(ss) and isComplete(ss, 0, K);
        };

        fmt::print ("{} ham/edit {}/{} → {}/{} ({}, {})\n", K, nc(ss), nce(ss), nc(ess), nce(ess), valid(ss), valid(ess));
    }
}
