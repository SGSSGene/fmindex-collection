#include <fmt/format.h>
#include <search_schemes/generator/pigeon.h>
#include <search_schemes/nodeCount.h>
#include <search_schemes/expectedNodeCount.h>
#include <search_schemes/expand.h>

int main(int argc, char** argv) {
    if (argc == 1) return 0;

    auto K = std::stod(argv[1]);
    {
        auto oss = search_schemes::generator::pigeon_opt(0, K);
        for (auto s : oss) {
            fmt::print("{}\n", fmt::join(s.pi, " "));
            fmt::print("{}\n", fmt::join(s.l, " "));
            fmt::print("{}\n\n", fmt::join(s.u, " "));
        }

        auto ss = search_schemes::expand(oss, 100);
        fmt::print("node count ham: {}\n",  search_schemes::nodeCount(ss, 4));
        fmt::print("node count edit: {}\n", search_schemes::nodeCountEdit(ss, 4));

        fmt::print("node count ham: {}\n",  search_schemes::expectedNodeCount(ss, 4, 3'000'000'000));
        fmt::print("node count edit: {}\n", search_schemes::expectedNodeCountEdit(ss, 4, 3'000'000'000));


        auto counts = search_schemes::expandCount(oss[0].pi.size(), 100);
        fmt::print("counts {}\n", fmt::join(counts, " "));

        for (auto s : ss) {
            fmt::print("{}\n", fmt::join(s.pi, " "));
            fmt::print("{}\n", fmt::join(s.l, " "));
            fmt::print("{}\n\n", fmt::join(s.u, " "));
        }
        auto ess = search_schemes::expandDynamic(oss, 100, 4, 3'000'000'000);
        for (auto s : ess) {
            fmt::print("{}\n", fmt::join(s.pi, " "));
            fmt::print("{}\n", fmt::join(s.l, " "));
            fmt::print("{}\n\n", fmt::join(s.u, " "));
        }
        fmt::print("node count ham: {}\n",  search_schemes::expectedNodeCount(ess, 4, 3'000'000'000));
        fmt::print("node count edit: {}\n", search_schemes::expectedNodeCountEdit(ess, 4, 3'000'000'000));


    }
}
