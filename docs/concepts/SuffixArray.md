<!--
    SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
    SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
    SPDX-License-Identifier: CC-BY-4.0
-->
# Suffix Array
## Suffix Arrays
Classes fulfilling this concept provide the member functions:

- `#!c++ auto value(size_t idx) const -> std::optional<std::tuple<uint64_t, uint64_t>>`

## Implementation

- `#!c++ fmindex_collection::CSA`
- `#!c++ fmindex_collection::DenseCSA`
