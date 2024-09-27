#pragma once

#include <format>
#include <iostream>
#include <string>
#include <vector>

struct BenchSize {
    struct Entry {
        std::string name;
        size_t      size;
        size_t      text_size;
        double      bits_per_char;
        double      relative{};
    };

    void addEntry(Entry e) {
        entries.push_back(e);
        entries.back().relative = double(entries.back().size) / entries[0].size * 100.;
    }

    std::vector<Entry> entries;

    ~BenchSize() {
        if (entries.empty()) return;
        auto lines = std::vector<std::string>{};
        lines.resize(entries.size());

        auto addColumn = [&]<typename T>(std::format_string<T&, size_t&> fmt, T BenchSize::Entry::*ptr) {
            size_t longestLine{};
            for (size_t i{0}; i < lines.size(); ++i) {
                size_t x = 1;
                auto t = std::format(fmt, entries[i].*ptr, x);
                longestLine = std::max(longestLine, t.size());
            }
            for (size_t i{0}; i < lines.size(); ++i) {
                auto t = std::format(fmt, entries[i].*ptr, longestLine);
                lines[i] = std::format("{} | {}", lines[i], t);
            }
        };

        addColumn("{:> {}.1f}%", &BenchSize::Entry::relative);
        addColumn("{:> {}}", &BenchSize::Entry::size);
        addColumn("{:> {}.3f}", &BenchSize::Entry::bits_per_char);
        addColumn("{:<{}}", &BenchSize::Entry::name);

        for (size_t i{0}; i < lines.size(); ++i) {
            std::cout << std::format("{}\n", lines[i]);
        }
    }
};
