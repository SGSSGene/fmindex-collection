// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: CC0-1.0

#include <fmindex-collection/fmindex-collection.h>
#include <fmt/format.h>

int main() {
    // your database/the data you want to search through, (value 0 should be avoided)
    auto reference = std::vector<std::vector<uint8_t>> {
        {1, 1, 1, 2, 2, 2, 3, 2, 4, 1, 1, 1}, // seqId 0
        {1, 2, 1, 2, 3, 4, 3},                // seqId 1
    };

    // Creating an bidirectional FM-Inddex over an alphabet with numbers 0-4
    auto index = fmc::BiFMIndex<5>{reference, /*samplingRate*/16, /*threadNbr*/1};

    // The stuff you are searching for
    auto query = std::vector<uint8_t>{2, 3};

    fmc::search</*.EditDistance=*/true>(index, std::span{&query, 1}, /*.errors=*/0, [&](size_t /*qidx*/, auto cursor, size_t errors) {
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
