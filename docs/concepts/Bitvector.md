<!--
    SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
    SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
    SPDX-License-Identifier: CC-BY-4.0
-->
# Bitvector

## Concept
The concept `fmindex_collection::Bitvector_c` models
a bitvector with rank support. Classes fulfilling this concept are guaranteed to provide following
functionality:

- class is `std::default_initializable<T>`
- class is `std::movable<T>`
- `#!c++ T(size_t length, Proj proj)`
- `#!c++ auto size() const -> size_t`
- `#!c++ auto symbol(size_t idx) const -> bool`
- `#!c++ auto rank(size_t idx) const -> bool`

## Examples usages
=== "Construction and Usage"
    Creating a bitvector of length 100,
    in which every second bit is set to one.
    ```c++
    auto bitvector = T{100, [](size_t idx) {
        return idx % 2 == 0;
    });

    // prints the number of bits
    std::cout << bitvector.size() << '\n';

    // prints value of the 51th bit
    std::cout << bitvector.symbol(50) << '\n';

    // prints how many ones there are in the first 50bit
    std::cout << bitvector.rank(50) << '\n';
    ```

=== "Construction V2"
    Creating a bitvector of length 100, where
    bits 0-49 are set to 1 and then move it around.
    ```c++
    auto bitvector = T{100, [](size_t idx) {
        return idx < 50;
    });

    T bitvector2 = std::move(bitvector);
    ```

## Implementations

- `#!c++ fmindex_collection::bitvector::Bitvector`
- `#!c++ fmindex_collection::bitvector::CompactBitvector`
- `#!c++ fmindex_collection::bitvector::CompactBitvector4Blocks`

## Stats
|                     Name | Space (bits per bit) | Description |
|:------------------------ | --------------------:| ----------- |
|`Bitvector`               |                1.375 | Consist of 3 independent arrays. 64bit values for bits, 8bit for blocks and 64bit for superblocks. Every 64bits a new block is started. Every 256bits a new superblock starts|
|`CompactBitvector`        |                1.333 | Consist of 3 interleaved arrays. Each DataBlock is 512 bits, and covers 6·64bits, each block is 9 bit and starts every 64bit, superblocks start ever 352 bits, 10bits are wasted for padding.|
|`CompactBitvector4Blocks` |                1.375 | Consist of 3 interleaved arrays. Each DataBlock is 352 bits, and covers 4·64bits, each block is 8 bit and starts every 64bit, superblocks start every 256 bits. |

![bit vector memory layout](Bitvector.png){ align=right }
