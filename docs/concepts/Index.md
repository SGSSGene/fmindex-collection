<!--
    SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
    SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
    SPDX-License-Identifier: CC-BY-4.0
-->
# FM-Index
We have to main indices. The `FMindex` and the `BiFMIndex`.

As a first parameter, they take a string implementation with an alphabet size.

## Example

Creates a BiFMIndex that covers the 4 base nucleotide A, C, G, T represented as 1, 2, 3, 4.
The BiFMIndex and the FMIndex require the `0` for special internal purposes.
It is possible to use it in a reference text, but it will not be considered during a search.


```c++
int main() {
    auto ref = std::vector<std::vector<uint8_t>>{
        {1, 1, 1, 2, 2, 2, 3, 2, 4, 1, 1, 1}, // seqId 0
        {1, 2, 1, 2, 3, 4, 3},                // seqId 1
    };

    // Pick String with Ranksupport implementation
    using String = fmc::string::FlattenedBitvectors_512_64k<5>; // 5 number of different characters

    // Creating an bidirectional FM-Inddex
    auto index = fmc::BiFMIndex<String>{reference, /*.samplingRate=*/16, /*.threadNbr=*/1};
}
```

# Advanced
Different FM-Indices classes exist:

- `#!cpp fmindex_collection::FMIndex<OccTable Table, typename TCSA>`
- `#!cpp fmindex_collection::BiFMIndex<OccTable Table, typename TCSA>`
- `#!cpp fmindex_collection::ReverseFMIndex<OccTable Table, typename TCSA>`
- `#!cpp fmindex_collection::RBiFMIndex<OccTable Table, typename TCSA>` (what does this one do?)
