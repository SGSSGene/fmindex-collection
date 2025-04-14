//SPDX-FileCopyrightText: 2024 Simon Gene Gottlieb
//SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <fmt/format.h>
#include <iostream>
#include <string>
#include <vector>

struct BenchSize {
    // # 4 columns:
    // - relative
    // - total size (in bits)
    // - bits per entry
    // - name
    std::vector<std::array<std::string, 5>> entries {
        {"relative", "size", "bits/bit", "overhead %", "name"}
    };

    double baseSize{0};

    struct Entry {
        std::string name;
        size_t      size;
        size_t      text_size;
        double      bits_per_char;
        double      relative{};
    };
    size_t firstEntrySize{};

    void addEntry(Entry e) {
        if (entries.size() == 1) {
            firstEntrySize = e.size;
        }

        double overheadInPercent = (e.bits_per_char - baseSize) / baseSize * 100.;

        entries.push_back({
            fmt::format("{:.1f}%", double(e.size) / firstEntrySize*100.),
            fmt::format(std::locale("en_US.UTF-8"), "{:L}", e.size),
            fmt::format("{:.3f}", e.bits_per_char),
            fmt::format("{:.3f}%", overheadInPercent),
            fmt::format("{}", e.name)
        });
    }

    ~BenchSize() {
        if (entries.size() < 2) return;
        fmt::print("\n");

        auto sizesPerColumn = std::array<size_t, 5>{};

        auto& c = sizesPerColumn;
        for (auto const& e : entries) {
            for (size_t i{0}; i < c.size(); ++i) {
                c[i] = std::max(c[i], e[i].size());
            }
        }

        for (size_t i{0}; i < entries.size(); ++i) {
            auto const& e = entries[i];
            auto const& sc = sizesPerColumn;
            if (i == 1) {
                fmt::print("|-{0:->{1}}-|-{0:->{2}}-|-{0:->{3}}-|-{0:-<{4}}-|-{0:-<{5}}\n", "", sc[0], sc[1], sc[2], sc[3], sc[4]);
            }

            fmt::print("| {: >{}} | {: >{}} | {: >{}} | {: >{}} | {: <{}} |\n", e[0], sc[0], e[1], sc[1], e[2], sc[2], e[3], sc[3], e[4], sc[4]);
        }
    }
};
