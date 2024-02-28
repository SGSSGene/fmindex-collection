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

| Functionality                                 | Description |
|-----------------------------------------------|-------------|
| class is `std::default_initializable<T>`      |             |
| class is `std::movable<T>`                    |             |
| `#!c++ T(size_t length, Proj proj)`           |             |
| `#!c++ auto size() const -> size_t`           | Returns the number of bits of this vector. |
| `#!c++ auto symbol(size_t idx) const -> bool` | Returns the bit at position `idx`. Parameter must be 0 ≤ idx < size(). |
| `#!c++ auto rank(size_t idx) const -> bool`   | Returns the number of ones in the first `idx` bits. Paremeter must be 0 ≤ idx ≤ size(). |

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

## Implementations - Description

- `#!c++ fmindex_collection::bitvector::Bitvector`

    Consist of 3 independent arrays. 64bit values for bits, 8bit for blocks and 64bit for superblocks. Every 64bits a new block is started. Every 256bits a new superblock starts.

- `#!c++ fmindex_collection::bitvector::CompactBitvector`

    Consist of 3 interleaved arrays. Each DataBlock is 512 bits, and covers 6·64bits, each block is 9 bit and starts every 64bit, superblocks start ever 352 bits, 10bits are wasted for padding.

- `#!c++ fmindex_collection::bitvector::CompactBitvector4Blocks`

    Consist of 3 interleaved arrays. Each DataBlock is 352 bits, and covers 4·64bits, each block is 8 bit and starts every 64bit, superblocks start every 256 bits.

- `#!c++ fmindex_collection::bitvector::SparseBLEBitvector`

    Sparse Block Length Encoded Bitvector. !TODO add paper reference

## Statistics
### Memory
Percentage refer to how many ones where present in the input data.

<div markdown class="compact_data_table">

|                         Name |   50%    |   25%    |   10%    |    5%    |  0.5%    |
|:-----------------------------|---------:|---------:|---------:|---------:|---------:|
|`Bitvector`                   |   1.38   |   1.38   |   1.38   |   1.38   |   1.38   |
|`CompactBitvector`            | **1.33** |   1.33   |   1.33   |   1.33   |   1.33   |
|`CompactBitvector4Blocks`     |   1.50   |   1.50   |   1.50   |   1.50   |   1.50   |
|`SparseBLEBitvector`  2       |   1.72   | **1.29** |   0.95   |   0.82   |   0.70   |
|`SparseBLEBitvector`  4       |   1.63   | **1.28** | **0.82** |   0.60   |   0.37   |
|`SparseBLEBitvector`  8       |   1.54   |   1.41   |   0.96   |   0.63   |   0.23   |
|`SparseBLEBitvector` 16       |   1.46   |   1.45   |   1.21   |   0.86   |   0.19   |
|`SparseBLEBitvector` 32       |   1.42   |   1.42   |   1.37   |   1.15   |   0.25   |
|`SparseBLEBitvector` 2/2      |   2.02   |   1.42   |   0.84   |   0.61   |   0.37   |
|`SparseBLEBitvector` 4/2      |   1.80   |   1.42   |   0.84   |   0.54   |   0.21   |
|`SparseBLEBitvector` 4/4      |   1.72   |   1.37   |   0.84   | **0.53** |   0.14   |
|`SparseBLEBitvector` 8/2      |   1.63   |   1.49   |   1.01   |   0.65   |   0.15   |
|`SparseBLEBitvector` 8/4      |   1.58   |   1.45   |   0.99   |   0.64   | **0.12** |
|`SparseBLEBitvector` 8/8      |   1.56   |   1.43   |   0.98   |   0.65   | **0.12** |
|`SparseBLEBitvector` 4/-2     |   2.02   |   1.42   |   0.84   |   0.61   |   0.37   |
|`SparseBLEBitvector` 8/-2     |   1.89   |   1.39   | **0.82** | **0.54** |   0.21   |
|`SparseBLEBitvector` 8/-4     |   1.80   |   1.42   |   0.84   | **0.54** |   0.21   |


</div>


### Run-time
Input data was 100'000'000 bits with 50% of them ones.

<div markdown class="compact_data_table">

|                         Name | Construct   |   symbol()   |   rank()    |
|:---------------------------- | -----------:|  ---------  :|  -------   :|
|`Bitvector`                   |  **1.97ns** |  **15.64ns** | **40.93ns** |
|`CompactBitvector`            |    1.98ns   |    40.51ns   |   90.74ns   |
|`CompactBitvector4Blocks`     |    2.39ns   |    28.72ns   |   47.80ns   |
|`SparseBLEBitvector` 4        |    4.14ns   |    27.86ns   |  183.78ns   |
|`SparseBLEBitvector` !4       |    3.58ns   |    40.56ns   |   56.90ns   |
|`SparseBLEBitvector` 8/4      |    4.01ns   |    45.92ns   |   78.69ns   |
|`SparseBLEBitvector` 8/-4     |    5.18ns   |    43.57ns   |  117.80ns   |

</div>


![bit vector memory layout](Bitvector.png){ align=right }
