// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: CC0-1.0

#include "../string/utils.h"
#include <fmindex-collection/fmindex/FMIndex.h>

TEST_CASE("benchmark fmindex on c'tor operation - 5 alphabet", "[fmindex][!benchmark][5][time][ctor]") {
    auto const& text = generateTexts<1, 4>();

    SECTION("benchmarking") {
        auto bench = ankerl::nanobench::Bench{};
        bench.title("c'tor()")
             .relative(true)
             .batch(text.size());

        call_with_templates<
            ALLSTRINGSWITHRANK(5)>([&]<typename String>() {
            if constexpr (std::same_as<String, fmindex_collection::string::Naive<5>>) {
                return;
            }

            auto string_name = getName<String>();
            INFO(string_name);

            bench.run(string_name, [&]() {
                using Index = fmindex_collection::FMIndex<String>;
                auto index = Index{text, /*.samplingRate=*/16, /*.threadNbr=*/1};
                ankerl::nanobench::doNotOptimizeAway(const_cast<Index const&>(index));
            });
        });
    }
}
