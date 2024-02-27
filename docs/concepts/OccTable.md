<!--
    SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
    SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
    SPDX-License-Identifier: CC-BY-4.0
-->
# Occurrence Table
These are mainly a wrapper around RankVector giving them additional functionality. More specific for the
FM-Index they provide the 'C' array. `GenericOccTable` is the direct translation of RankVector and adding a 'C' array.
All others are special implementations that haven't been refactored yet.


## OccTable (`fmindex_collection::OccTable`)
Classes fulfilling this concept provide the member functions:

- `#!c++ auto rank(size_t idx, uint8_t symb) const -> size_t`
- `#!c++ auto prefix_rank(size_t idx, uint8_t symb) const -> size_t`
- `#!c++ auto all_ranks(size_t) const -> std::tuple<std::array<size_t, T::Sigma>, std::array<size_t, T::Sigma>>`
- `#!c++ auto size() const -> size_t`

## Implementation

| occ table (`fmindex_collection::`)         | Description |
|:-------------------------------------------|--------------|
| `#!c++ GenericOccTable<SymbolVector Vector>`                            | Wrapper around a rank supported vector |
| `#!c++ occtable::compactBitvectorPrefix::OccTable<size_t TSigma>`       | uses internal bit vectors representing prefix values |
| `#!c++ occtable::eprV8::OccTable<size_t TSigma>`                        | Variable sized alphabet vector (`TSigma` has no functionality) |
| `#!c++ occtable::interleavedPrefix::OccTablee<size_t TSigma>`           | uses internal bit vectors representing prefix values |
