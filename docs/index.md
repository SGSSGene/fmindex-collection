<!--
    SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
    SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
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
### Quick start
```c++
#include <fmindex-collection/fmindex-collection.h>
#include <fmt/format.h>

namespace fmc = fmindex_collection;

int main() {
    // your database/the data you want to search through, (value 0 should be avoided)
    auto reference = std::vector<std::vector<uint8_t>> {
        {1, 1, 1, 2, 2, 2, 3, 2, 4, 1, 1, 1}, // seqId 0
        {1, 2, 1, 2, 3, 4, 3},                // seqId 1
    };

    // Pick String with Ranksupport implementation
    using String = fmc::string::FlattenedBitvectors_512_64k<5>; // 5 number of different characters

    // Creating an bidirectional FM-Inddex
    auto index = fmc::BiFMIndex<String>{reference, /*samplingRate*/16, /*threadNbr*/1};

    // The stuff you are searching for
    auto query = std::vector<uint8_t>{2, 3};

    fmc::search</*EditDistance=*/true>(index, query, /*.errors=*/0, [&](auto cursor, size_t errors) {
        fmt::print("found {} results with {} errors\n", cursor.count(), errors);
        for (auto i : cursor) {
            // index.locate(i) can only find positions hit by the sampling rate. How many position this is off, is indicated by the offset value
            auto [entry, offset] = index.locate(i);
            auto [seqId, pos] = entry; // tuple of the sequence id and the position inside the sequence
            fmt::print("seqId/pos: {}/{}\n", seqId, pos+offset);
        }
    });
    return 0;
}
```
