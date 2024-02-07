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
- `#!c++ fmindex_collection::rankvector::CompactBitvector` (obsolete?)
- `#!c++ fmindex_collection::rankvector::DenseEPRV6<size_t TSigma>`
- `#!c++ fmindex_collection::rankvector::EPRV3_8<size_t TSigma>`
- `#!c++ fmindex_collection::rankvector::EPRV3_16<size_t TSigma>`
- `#!c++ fmindex_collection::rankvector::EPRV3_32<size_t TSigma>`
- `#!c++ fmindex_collection::rankvector::EPRV4<size_t TSigma>`
- `#!c++ fmindex_collection::rankvector::EPRV5<size_t TSigma>`
- `#!c++ fmindex_collection::rankvector::InterleavedBitvector8<size_t TSigma>`
- `#!c++ fmindex_collection::rankvector::InterleavedBitvector16<size_t TSigma>`
- `#!c++ fmindex_collection::rankvector::InterleavedBitvector32<size_t TSigma>`
- `#!c++ fmindex_collection::rankvector::InterleavedBitvector8Aligned<size_t TSigma>`
- `#!c++ fmindex_collection::rankvector::InterleavedBitvector16Aligned<size_t TSigma>`
- `#!c++ fmindex_collection::rankvector::InterleavedBitvector32Aligned<size_t TSigma>`
- `#!c++ fmindex_collection::rankvector::InterleavedEPR8<size_t TSigma>`
- `#!c++ fmindex_collection::rankvector::InterleavedEPR16<size_t TSigma>`
- `#!c++ fmindex_collection::rankvector::InterleavedEPR32<size_t TSigma>`
- `#!c++ fmindex_collection::rankvector::InterleavedEPR8Aligned<size_t TSigma>`
- `#!c++ fmindex_collection::rankvector::InterleavedEPR16Aligned<size_t TSigma>`
- `#!c++ fmindex_collection::rankvector::InterleavedEPR32Aligned<size_t TSigma>`
- `#!c++ fmindex_collection::rankvector::InterleavedEPRV2_8<size_t TSigma>`
- `#!c++ fmindex_collection::rankvector::InterleavedEPRV2_16<size_t TSigma>`
- `#!c++ fmindex_collection::rankvector::InterleavedEPRV2_32<size_t TSigma>`
- `#!c++ fmindex_collection::rankvector::InterleavedEPRV2_8Aligned<size_t TSigma>`
- `#!c++ fmindex_collection::rankvector::InterleavedEPRV2_16Aligned<size_t TSigma>`
- `#!c++ fmindex_collection::rankvector::InterleavedEPRV2_32Aligned<size_t TSigma>`
- `#!c++ fmindex_collection::rankvector::InterleavedWavelet<size_t TSigma>`
- `#!c++ fmindex_collection::rankvector::MultiBitvector<size_t TSigma>`
- `#!c++ fmindex_collection::rankvector::Naive<size_t TSigma>`
- `#!c++ fmindex_collection::rankvector::RLE<size_t TSigma, size_t encodingBlockSize>` (some parameters are missing)
- `#!c++ fmindex_collection::rankvector::rRLE<size_t TSigma, size_t encodingBlockSize>` (some parameters are missing)
- `#!c++ fmindex_collection::rankvector::Sdsl_wt_bldc<size_t TSigma>`
- `#!c++ fmindex_collection::rankvector::Sdsl_wt_epr<size_t TSigma>` (something is broken here)
- `#!c++ fmindex_collection::rankvector::Wavelet<size_t TSigma>`


## Stats
|                                           Name | Space (bits per bit) |
|:---------------------------------------------- | --------------------:|
- `CompactBitvector`                             |                !TODO |
- `DenseEPRV6<size_t TSigma>`                    |                !TODO |
- `EPRV3_8<size_t TSigma>`                       |                !TODO |
- `EPRV3_16<size_t TSigma>`                      |                !TODO |
- `EPRV3_32<size_t TSigma>`                      |                !TODO |
- `EPRV4<size_t TSigma>`                         |                !TODO |
- `EPRV5<size_t TSigma>`                         |                !TODO |
- `InterleavedBitvector8<size_t TSigma>`         |              n·log n |
- `InterleavedBitvector16<size_t TSigma>`        |                !TODO |
- `InterleavedBitvector32<size_t TSigma>`        |                !TODO |
- `InterleavedBitvector8Aligned<size_t TSigma>`  |                !TODO |
- `InterleavedBitvector16Aligned<size_t TSigma>` |                !TODO |
- `InterleavedBitvector32Aligned<size_t TSigma>` |                !TODO |
- `InterleavedEPR8<size_t TSigma>`               |                !TODO |
- `InterleavedEPR16<size_t TSigma>`              |                !TODO |
- `InterleavedEPR32<size_t TSigma>`              |                !TODO |
- `InterleavedEPR8Aligned<size_t TSigma>`        |                !TODO |
- `InterleavedEPR16Aligned<size_t TSigma>`       |                !TODO |
- `InterleavedEPR32Aligned<size_t TSigma>`       |                !TODO |
- `InterleavedEPRV2_8<size_t TSigma>`            |                !TODO |
- `InterleavedEPRV2_16<size_t TSigma>`           |                !TODO |
- `InterleavedEPRV2_32<size_t TSigma>`           |                !TODO |
- `InterleavedEPRV2_8Aligned<size_t TSigma>`     |                !TODO |
- `InterleavedEPRV2_16Aligned<size_t TSigma>`    |                !TODO |
- `InterleavedEPRV2_32Aligned<size_t TSigma>`    |                !TODO |
- `InterleavedWavelet<size_t TSigma>`            |                !TODO |
- `MultiBitvector<size_t TSigma>`                |            σ·n·1.375 |
- `Naive<size_t TSigma>`                         |            σ·n·log n |
- `RLE<size_t TSigma, size_t encodingBlockSize>` |                !TODO |
- `rRLE<size_t TSigma, size_t encodingBlockSize>`|                !TODO |
- `Sdsl_wt_bldc<size_t TSigma>`                  |                !TODO |
- `Sdsl_wt_epr<size_t TSigma>`                   |                !TODO |
- `Wavelet<size_t TSigma>`                       |                !TODO |
