// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#include <cstdio>
#include <fmindex-collection/search_scheme/generator/all.h>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <unordered_set>

namespace ss = fmc::search_scheme;

void help() {
    fmt::print("Usage:\n"
                "./search_scheme_generator <generator> <k>\n\n"
                "generators:\n");

    for (auto const& [key, value] : ss::generator::all) {
        fmt::print("- {}\n", key);
    }
}

int main(int argc, char const* const* argv) {

    if (argc < 3 || std::string_view{argv[1]} == "--help") {
        help();
        return 0;
    }
    try {
        auto generator_name = std::string{argv[1]};
        auto k              = std::stoi(std::string{argv[2]});

        auto iter = ss::generator::all.find(generator_name);
        if (iter == ss::generator::all.end()) {
            throw std::runtime_error("unknown search scheme generetaror \"" + generator_name + "\"");
        }

        auto oss = iter->second.generator(0, k, 0, 0); //!TODO last two parameters of second are not being used
        fmt::print("\\item \\(k={}\\): ", k);
        fmt::print("\\{{");
        for (size_t i{0}; i < oss.size(); ++i) {
            auto const& search = oss[i];
            if (i > 0) fmt::print(", ");
            fmt::print("\\(({}, {}, {})\\)", fmt::join(search.pi, ""), fmt::join(search.l, ""), fmt::join(search.u, ""));
        }
        fmt::print("\\}}\n");
//        auto ess = ss::expand(oss, len);


    } catch(std::exception const& e) {
        fmt::print("{}\n===\n\n", e.what());
        help();
    }
    return 0;
}
