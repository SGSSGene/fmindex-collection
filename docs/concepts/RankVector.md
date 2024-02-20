<!--
    SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
    SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
    SPDX-License-Identifier: CC-BY-4.0
-->
# RankVector

## Concept
The concept `fmindex_collection::RankVector` models
a vector that stores elements of type `uint8_t` and provides a `rank` function for each symbol.

Classes fulfilling this concept are guaranteed to provide following
functionality:

- class is `#!c++ std::default_initializable<T>`
- class is `#!c++ std::movable<T>`
- `#!c++ constexpr size_t Sigma`
- `#!c++ T{std::span<uint8_t> data)`
- `#!c++ auto size() const -> size_t`
- `#!c++ auto rank(size_t idx, uint8_t symb) const -> size_t`
- `#!c++ auto prefix_rank(size_t idx, uint8_t symb) const -> size_t`
- `#!c++ auto all_ranks(size_t idx) const -> std::array<size_t, T::Sigma>`
- `#!c++ auto all_ranks_and_prefix_ranks(size_t idx) const -> std::tuple<std::array<size_t, T::Sigma>, std::array<size_t, T::Sigma>>`

## Examples usages
=== "Construction and Usage"
    Creating a vector with rank support of length 100.
    ```c++
    auto data = std::vector<uint8_t>{3, 2, 3, 2, 1, 1, 0, 3, 4, 5};
    auto vector = T{data};

    // prints the number of bits
    std::cout << vector.size() << '\n';

    // prints value of 5th bit (prints '1')
    std::cout << vector.symbol(4) << '\n';

    // prints how number of '2' in the first 4 numbers (prints '2')
    std::cout << bitvector.rank(4, 2) << '\n';
    ```

## Implementation

All implementation are located inside the namespace `fmindex_collection::rankvector`.

<div markdown class="rankvector_table">

| rank vector                                               | σ=4 | σ=5 | σ=6 | σ=21 | σ=256(!TODO for 1GB of data, currently only 100MB)| Description |
|:----------------------------------------------------------|------:|---:|---:|---:|----:|--------------|
| `#!c++ Naive`                              |    256| 320| 384|1344|16384| Stores every thing in `std::vector<size_t>`, requires O(n·log n). |
| `#!c++ MultiBitvector`                     |    5.5| 6.9| 8.3|28.9|352.0| Standard FM-Index implementation. Using a bitvector for every alphabet character. |
| `#!c++ InterleavedBitvector8`              |    5.5| 6.9| 8.3|28.9|352.0| Interleaving blocks and bits. Using 8 bits for storing block values. New super-block every 256bits.     |
| `#!c++ InterleavedBitvector16`             |    5.0| 6.3| 7.5|26.3|320.3| Interleaving blocks and bits. Using 16 bits for storing block values. New super-block every 65536 bits. |
| `#!c++ InterleavedBitvector32`             |    6.0| 7.5| 9.0|31.5|384.0| Interleaving blocks and bits. Using 32 bits for storing block values. New super-block every 2^32 bits.  |
| `#!c++ InterleavedBitvector8Aligned`       |       |    |    |    |     | Same as `InterleavedBitvector` but the bits and blocks sections are aligned to 64bit. |
| `#!c++ InterleavedBitvector16Aligned`      |       |    |    |    |     | Same as `InterleavedBitvector` but the bits and blocks sections are aligned to 64bit. |
| `#!c++ InterleavedBitvector32Aligned`      |       |    |    |    |     | Same as `InterleavedBitvector` but the bits and blocks sections are aligned to 64bit. |
| `#!c++ InterleavedEPR8`                    |    4.0| 6.2| 6.9|24.7|328.0| Similar to `InterleavedBitvector, but bits are encoded according to their BWT representation. |
| `#!c++ InterleavedEPR16`                   |    4.0| 6.9| 7.6|33.4|520.3| See `InterleavedEPR8`. |
| `#!c++ InterleavedEPR32`                   |    6.0|10.7|12.2|61.3| 1032| See `InterleavedEPR8`. |
| `#!c++ InterleavedEPR8Aligned`             |       |    |    |    |     | See `InterleavedEPR8`. |
| `#!c++ InterleavedEPR16Aligned`            |       |    |    |    |     | See `InterleavedEPR8`. |
| `#!c++ InterleavedEPR32Aligned`            |       |    |    |    |     | See `InterleavedEPR8`. |
| `#!c++ InterleavedEPRV2_8`                 |    4.0| 5.3| 5.5|13.3|104.0| Similar to `InterleavedBitvector8`, but custom bit shuffling. Faster and better use of space compared to EPR. |
| `#!c++ InterleavedEPRV2_16`                |    3.0| 5.0| 5.0|11.0| 72.3| See `InterleavedEPRV2_8`. |
| `#!c++ InterleavedEPRV2_32`                |    4.0| 6.0| 6.0|16.0|136.0| See `InterleavedEPRV2_8`. |
| `#!c++ InterleavedEPRV2_8Aligned`          |       |    |    |    |     | See `InterleavedEPRV2_8`. |
| `#!c++ InterleavedEPRV2_16Aligned`         |       |    |    |    |     | See `InterleavedEPRV2_8`. |
| `#!c++ InterleavedEPRV2_32Aligned`         |       |    |    |    |     | See `InterleavedEPRV2_8`. |
| `#!c++ EPRV3_8`                            |    3.5| 4.9| 5.3|12.9|104.0| Similar to `InterlavedEPRV2_8` but no interleaving. |
| `#!c++ EPRV3_16`                           |    3.0| 4.3| 4.5|10.3| 72.3| See `EPRV3_8` |
| `#!c++ EPRV3_32`                           |    4.0| 5.5| 6.0|15.5|136.0| See `EPRV3_8` |
| `#!c++ EPRV4`                              |    2.8| 3.9| 4.1| 8.9| 56.1| Similar to  `EPRV3` but using two additional layers of blocks. (8bit, 16bit and 32bit). |
| `#!c++ EPRV5`                              |    2.8| 3.9| 4.1| 9.0| 56.3| Similar to `EPRV3` but using one additional layer of blocks. (8bit and 16bit). |
| `#!c++ DenseEPRV6`                         |    2.8| 3.9| 4.1| 8.9| 56.1| Same as `EPRV5` but using a dense vector for the super blocks instead of 64bit values. (Not worth it). |
| `#!c++ InterleavedEPRV7`                   |    2.8| 3.9| 4.1| 8.9| 56.3| Similar to `EPRV5` but blocks and lowest layer are interleaved. |
| `#!c++ InterleavedWavelet`                 |    5.0| 5.5| 6.0|15.5|137.0| Special version of Wavelet-trees (!TODO I think this was a combination fo EPR and Wavelet, not sure)|
| `#!c++ RLE`                                |       |    |    |    |     | (!TODO some parameters are missing) |
| `#!c++ rRLE`                               |       |    |    |    |     | (!TODO some parameters are missing) |
| `#!c++ Sdsl_wt_bldc`                       |    2.5| 3.0| 3.3| 5.6| 10.0| Wrapper of the Wavelet-Tree implementation of the SDSL library. |
| `#!c++ Sdsl_wt_epr`                        |       |    |    |    |     | Wrapper of the EPR implementation of the SDSL library. (!TODO something is broken) |
| `#!c++ Wavelet`                            |    4.0| 4.0| 4.0| 6.7| 12.0| Custom Wavelet-Tree implementation. |
| `#!c++ CompactBitvector`                   |       |    |    |    |     | Special case for Sigma 2, using `bitvector::CompactBitvector` as implementation. |
</div>
