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
            ALLSTRINGSWITHRANK>([&]<template <size_t> typename String>() {
            auto string_name = getName<String<5>>();
            INFO(string_name);

            bench.run(string_name, [&]() {
                using Index = fmindex_collection::FMIndex<5, String>;
                auto index = Index{text, /*.samplingRate=*/16, /*.threadNbr=*/1};
                ankerl::nanobench::doNotOptimizeAway(const_cast<Index const&>(index));
            });
        });
    }
}
