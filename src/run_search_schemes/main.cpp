#include <fmt/format.h>
#include <search_schemes/generator/pigeon.h>
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
    }
}
