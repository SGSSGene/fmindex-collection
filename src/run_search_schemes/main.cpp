#include <fmt/format.h>
#include <search_schemes/generator/pigeon.h>
#include <search_schemes/nodeCount.h>
#include <search_schemes/expectedNodeCount.h>
#include <search_schemes/expand.h>

int main(int argc, char** argv) {
    if (argc == 1) return 0;

    auto K = std::stod(argv[1]);
    {
        auto ss = search_schemes::generator::pigeon_opt(0, K);
        for (auto s : ss) {
            fmt::print("{}\n", fmt::join(s.pi, ""));
            fmt::print("{}\n", fmt::join(s.l, ""));
            fmt::print("{}\n\n", fmt::join(s.u, ""));
        }

        ss = search_schemes::expand(ss, 100);
        fmt::print("node count ham: {}\n",  search_schemes::nodeCount(ss, 4));
        fmt::print("node count edit: {}\n", search_schemes::nodeCountEdit(ss, 4));

        fmt::print("node count ham: {}\n",  search_schemes::expectedNodeCount(ss, 4, 3'000'000'000));
        fmt::print("node count edit: {}\n", search_schemes::expectedNodeCountEdit(ss, 4, 3'000'000'000));

    }
}
