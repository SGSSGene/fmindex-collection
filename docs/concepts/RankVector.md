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

| rank vector (`fmindex_collection::rankvecto::`) | Description |
|:-------------------------------------------|--------------|
| `#!c++ Naive`                              | Stores every thing in `std::vector<size_t>`, requires O(n·log n). |
| `#!c++ MultiBitvector`                     | Standard FM-Index implementation. Using a bitvector for every alphabet character. |
| `#!c++ InterleavedBitvector8`              | Interleaving blocks and bits. Using 8 bits for storing block values. New super-block every 256bits.     |
| `#!c++ InterleavedBitvector16`             | Interleaving blocks and bits. Using 16 bits for storing block values. New super-block every 65536 bits. |
| `#!c++ InterleavedBitvector32`             | Interleaving blocks and bits. Using 32 bits for storing block values. New super-block every 2^32 bits.  |
| `#!c++ InterleavedBitvector8Aligned`       | Same as `InterleavedBitvector` but the bits and blocks sections are aligned to 64bit. |
| `#!c++ InterleavedBitvector16Aligned`      | Same as `InterleavedBitvector` but the bits and blocks sections are aligned to 64bit. |
| `#!c++ InterleavedBitvector32Aligned`      | Same as `InterleavedBitvector` but the bits and blocks sections are aligned to 64bit. |
| `#!c++ InterleavedEPR8`                    | Similar to `InterleavedBitvector, but bits are encoded according to their BWT representation. |
| `#!c++ InterleavedEPR16`                   | See `InterleavedEPR8`. |
| `#!c++ InterleavedEPR32`                   | See `InterleavedEPR8`. |
| `#!c++ InterleavedEPR8Aligned`             | See `InterleavedEPR8`. |
| `#!c++ InterleavedEPR16Aligned`            | See `InterleavedEPR8`. |
| `#!c++ InterleavedEPR32Aligned`            | See `InterleavedEPR8`. |
| `#!c++ InterleavedEPRV2_8`                 | Similar to `InterleavedBitvector8`, but custom bit shuffling. Faster and better use of space compared to EPR. |
| `#!c++ InterleavedEPRV2_16`                | See `InterleavedEPRV2_8`. |
| `#!c++ InterleavedEPRV2_32`                | See `InterleavedEPRV2_8`. |
| `#!c++ InterleavedEPRV2_8Aligned`          | See `InterleavedEPRV2_8`. |
| `#!c++ InterleavedEPRV2_16Aligned`         | See `InterleavedEPRV2_8`. |
| `#!c++ InterleavedEPRV2_32Aligned`         | See `InterleavedEPRV2_8`. |
| `#!c++ EPRV3_8`                            | Similar to `InterlavedEPRV2_8` but no interleaving. |
| `#!c++ EPRV3_16`                           | See `EPRV3_8` |
| `#!c++ EPRV3_32`                           | See `EPRV3_8` |
| `#!c++ EPRV4`                              | Similar to  `EPRV3` but using two additional layers of blocks. (8bit, 16bit and 32bit). |
| `#!c++ EPRV5`                              | Similar to `EPRV3` but using one additional layer of blocks. (8bit and 16bit). |
| `#!c++ DenseEPRV6`                         | Same as `EPRV5` but using a dense vector for the super blocks instead of 64bit values. (Not worth it). |
| `#!c++ InterleavedEPRV7`                   | Similar to `EPRV5` but blocks and lowest layer are interleaved. |
| `#!c++ InterleavedWavelet`                 | Special version of Wavelet-trees (!TODO I think this was a combination fo EPR and Wavelet, not sure)|
| `#!c++ RLE`                                | (!TODO some parameters are missing) |
| `#!c++ rRLE`                               | (!TODO some parameters are missing) |
| `#!c++ Sdsl_wt_bldc`                       | Wrapper of the Wavelet-Tree implementation of the SDSL library. |
| `#!c++ Sdsl_wt_epr`                        | Wrapper of the EPR implementation of the SDSL library. (!TODO something is broken) |
| `#!c++ Wavelet`                            | Custom Wavelet-Tree implementation. |
| `#!c++ CompactBitvector`                   | Special case for Sigma 2, using `bitvector::CompactBitvector` as implementation. |

</div>

## Statistics
### Memory

!TODO, the σ=256 ran only on 100MB of data, not 1GB of data.
<div markdown class="compact_data_table">

| rank vector                                |     σ=4 |    σ=5 |   σ=6 |   σ=21 |   σ=256 |
|:-------------------------------------------|--------:|-------:|------:|-------:|--------:|
| `#!c++ Naive`                              |    256  |   320  |  384  |  1344  |  16384  |
| `#!c++ MultiBitvector`                     |    5.5  |   6.9  |  8.3  |  28.9  |  352.0  |
| `#!c++ MultiBitvector` Sparse 2            |    5.2  |   5.9  |  6.7  |  17.1  |  178.8  |
| `#!c++ MultiBitvector` Sparse 4            |    5.1  |   5.8  |  6.3  |  12.3  |   93.5  |
| `#!c++ MultiBitvector` Sparse 8            |    5.6  |   6.6  |  7.4  |  12.9  |   54.9  |
| `#!c++ MultiBitvector` Sparse 16           |         |        |       |        | **43.4**|
| `#!c++ MultiBitvector` Sparse 32           |         |        |       |        |   52.5  |
| `#!c++ MultiBitvector` Sparse 64           |         |        |       |        |   83.5  |
| `#!c++ MultiBitvector` Sparse 2/2          |    6.0  |   7.2  |  8.3  |  24.2  |  266.8  |
| `#!c++ MultiBitvector` Sparse 4/2          |    5.2  |   6.0  |  6.8  |  15.7  |  137.5  |
| `#!c++ MultiBitvector` Sparse 4/4          |    5.2  |   6.0  |  6.7  |  14.1  |  115.5  |
| `#!c++ MultiBitvector` Sparse 4/-2         |    5.7  |   6.2  |  6.7  |  12.5  |   93.5  |
| `#!c++ MultiBitvector` Sparse 8/-2         |    5.6  |   6.2  |  6.7  |  11.0  |   52.2  |
| `#!c++ MultiBitvector` Sparse 8/-4         |    5.7  |   6.4  |  6.9  |  11.1  |   52.2  |
| `#!c++ MultiBitvector` Sparse 16/-4        |         |        |       |  10.8  | **32.9**|
| `#!c++ InterleavedBitvector8`              |    5.5  |   6.9  |  8.3  |  28.9  |  352.0  |
| `#!c++ InterleavedBitvector16`             |    5.0  |   6.3  |  7.5  |  26.3  |  320.3  |
| `#!c++ InterleavedBitvector32`             |    6.0  |   7.5  |  9.0  |  31.5  |  384.0  |
| `#!c++ InterleavedBitvector8Aligned`       |         |        |       |        |         |
| `#!c++ InterleavedBitvector16Aligned`      |         |        |       |        |         |
| `#!c++ InterleavedBitvector32Aligned`      |         |        |       |        |         |
| `#!c++ InterleavedEPR8`                    |    4.0  |   6.2  |  6.9  |  24.7  |  328.0  |
| `#!c++ InterleavedEPR16`                   |    4.0  |   6.9  |  7.6  |  33.4  |  520.3  |
| `#!c++ InterleavedEPR32`                   |    6.0  |  10.7  | 12.2  |  61.3  |   1032  |
| `#!c++ InterleavedEPR8Aligned`             |         |        |       |        |         |
| `#!c++ InterleavedEPR16Aligned`            |         |        |       |        |         |
| `#!c++ InterleavedEPR32Aligned`            |         |        |       |        |         |
| `#!c++ InterleavedEPRV2_8`                 |    4.0  |   5.3  |  5.5  |  13.3  |  104.0  |
| `#!c++ InterleavedEPRV2_16`                |    3.0  |   5.0  |  5.0  |  11.0  |   72.3  |
| `#!c++ InterleavedEPRV2_32`                |    4.0  |   6.0  |  6.0  |  16.0  |  136.0  |
| `#!c++ InterleavedEPRV2_8Aligned`          |         |        |       |        |         |
| `#!c++ InterleavedEPRV2_16Aligned`         |         |        |       |        |         |
| `#!c++ InterleavedEPRV2_32Aligned`         |         |        |       |        |         |
| `#!c++ EPRV3_8`                            |    3.5  |   4.9  |  5.3  |  12.9  |  104.0  |
| `#!c++ EPRV3_16`                           |    3.0  |   4.3  |  4.5  |  10.3  |   72.3  |
| `#!c++ EPRV3_32`                           |    4.0  |   5.5  |  6.0  |  15.5  |  136.0  |
| `#!c++ EPRV4`                              |    2.8  |   3.9  |  4.1  |   8.9  |   56.1  |
| `#!c++ EPRV5`                              |    2.8  |   3.9  |  4.1  |   9.0  |   56.3  |
| `#!c++ DenseEPRV6`                         |    2.8  |   3.9  |  4.1  |   8.9  |   56.1  |
| `#!c++ InterleavedEPRV7`                   |  **2.8**| **3.9**|**4.1**| **8.9**| **56.3**|
| `#!c++ InterleavedWavelet`                 |    5.0  |   5.5  |  6.0  |  15.5  |  137.0  |
| `#!c++ RLE`                                |         |        |       |        |         |
| `#!c++ rRLE`                               |         |        |       |        |         |
| `#!c++ Sdsl_wt_bldc`                       |  **2.5**| **3.0**|**3.3**| **5.6**| **10.0**|
| `#!c++ Sdsl_wt_epr`                        |         |        |       |        |         |
| `#!c++ Wavelet`                            |    4.0  |   4.0  |  4.0  |   6.7  |   12.0  |
| `#!c++ CompactBitvector`                   |         |        |       |        |         |

</div>

### Run-time
Input data was 100'000 characters in range of 1-4

<div markdown class="compact_data_table">

|                                                                 Name | Construct(ns)| symbol() (ns) |  rank() (ns)   |
|:-------------------------------------------------------------------- |-------------:|--------------:|---------------:|
| Naive<256>                                                           |   252.89ns  |        na        |          na    |
| MultiBitvector<256, bitvector::Bitvector>                            |   218.50ns  |     39.99      |       42.66  |
| MultiBitvector<256, bitvector::CompactBitvector>                     |   217.75ns  |     44.17      |       27.69  |
| MultiBitvector<256, bitvector::CompactBitvector4Blocks>              |   220.26ns  |     44.54      |       28.03  |
| MultiBitvector<256, bitvector::SparseBLEBitvector<>>  (defect)       |   453.69ns  |    132.20      |       62.22  |
| InterleavedBitvector8<256>                                           |     7.11ns  |     29.57      |       26.70  |
| InterleavedBitvector16<256>                                          |     5.15ns  |     28.40      |     **24.53**|
| InterleavedBitvector32<256>                                          |     6.62ns  |     28.52      |       25.12  |
| InterleavedBitvector8Aligned<256>                                    |     6.91ns  |     28.17      |       26.48  |
| InterleavedBitvector16Aligned<256>                                   |     5.07ns  |     28.03      |     **24.44**|
| InterleavedBitvector32Aligned<256>                                   |     6.46ns  |     28.39      |       25.06  |
| InterleavedEPR8<256>                                                 |     7.59ns  |   **14.83**    |       50.80  |
| InterleavedEPR16<256>                                                |    18.38ns  |   **14.86**    |       42.69  |
| InterleavedEPR32<256>                                                |    39.76ns  |  16158.44      |    19360.74  |
| InterleavedEPR8Aligned<256>                                          |    12.18ns  |     16.11      |       51.17  |
| InterleavedEPR16Aligned<256>                                         |    21.47ns  |     16.24      |       44.59  |
| InterleavedEPR32Aligned<256>                                         |    41.81ns  |  15377.29      |    18771.09  |
| InterleavedEPRV2_8<256>                                              |     3.42ns  |     39.15      |       40.15  |
| InterleavedEPRV2_16<256>                                             |   **2.55ns**|     38.89      |       43.58  |
| InterleavedEPRV2_32<256>                                             |     3.73ns  |     47.48      |       55.08  |
| InterleavedEPRV2_8Aligned<256>                                       |     3.43ns  |     38.04      |       44.63  |
| InterleavedEPRV2_16Aligned<256>                                      |   **2.43ns**|     37.40      |       43.29  |
| InterleavedEPRV2_32Aligned<256>                                      |     3.91ns  |     37.21      |       52.23  |
| EPRV3_8<256>                                                         |     3.65ns  |     36.76      |       61.26  |
| EPRV3_16<256>                                                        |     2.81ns  |     36.62      |       61.12  |
| EPRV3_32<256>                                                        |     2.92ns  |     36.67      |       46.78  |
| EPRV4<256>                                                           |     4.14ns  |     36.75      |       48.82  |
| EPRV5<256>                                                           |     4.07ns  |     36.28      |       46.65  |
| DenseEPRV6<256>                                                      |     3.97ns  |     37.59      |       65.85  |
| InterleavedEPRV7<256>                                                |     3.99ns  |     38.32      |       45.06  |
| InterleavedWavelet<256>                                              |    35.91ns  |    296.69      |      146.26  |
| Wavelet<256>                                                         |    27.85ns  |   1408.59      |      615.17  |
| RLE<256, 4, InterleavedBitvector16<256>, InterleavedBitvector16<256>>|     6.28ns  |     54.24      |       82.94  |
| Sdsl_wt_bldc<256>                                                    |    10.83ns  |     76.36      |      108.03  |

</div>
