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
|                     Name | Space (bits per bit) |
|:------------------------ | --------------------:|
|`Bitvector`               |                1.375 |
|`CompactBitvector`        |                1.333 |
|`CompactBitvector4Blocks` |                1.375 |
