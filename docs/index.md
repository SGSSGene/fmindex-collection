<!--
    SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
    SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
    SPDX-License-Identifier: CC-BY-4.0
-->
# Introduction

**FMIndex Collection** is written in modern C++20 and provides a set of concepts, classes and functions required
to provide FM-Indices and their internal support structures.

## Functionality
- fast and compact bitvectors with rank support
- fast and compact strings with rank support
- (bidirectional) FM-Indices
- search schemes generator
- several search algorithms (with edit distance or hamming distance)


## Dependencies
- Uses [libsais](https://github.com/IlyaGrebnov/libsais), Thank you [@IlyaGbrenov](https://github.com/IlyaGrebnov)

## Usage
### CPM
```cmake
CPMAddPackage(
  URI SGSSGene/fmindex_collection@0.0.0 # put newest version here
)
...
target_link_libraries(${PROJECT_NAME}
    fmindex_collection::fmindex_collection
)
```
### C++
```c++
#include <fmindex-collection/fmindex-collection.h>
#include <fmt/format.h>
void someFunction() {
    // your database/the data you want to search through
    auto reference = std::vector<std::vector<uint8_t>> {
        {1, 1, 1, 2, 2, 2, 3, 2, 4, 1, 1, 1},
        {1, 2, 1, 2, 3, 4, 3},
    };

    // Pick String with Ranksupport implementation
    using String = string::InterleavedBitvector16<5>; // largest character + 1 (not allowed to have 0)

    // Creating an FM-Inddex
    auto index = FMIndex<String>{reference, /*samplingRate*/16, /*threadNbr*/1};

    // The stuff you are searching for
    auto queries = std::vector<std::vector<uint8_t>>{{1, 3}, {4, 2}};

    search_backtracking::search(index, queries, /*.numberOfAllowedErrors=*/0, [&](size_t queryId, auto cursor, size_t errors) {
        (void)errors;
        fmt::print("found something {} {}\n", queryId, cursor.count());
        for (auto i : cursor) {
            auto [chr, pos] = index.locate(i);
            fmt::print("chr/pos: {}/{}\n", chr, pos);
        }
    });
```
