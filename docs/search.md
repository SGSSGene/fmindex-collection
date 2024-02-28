<!--
    SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
    SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
    SPDX-License-Identifier: CC-BY-4.0
-->
# Search

If you don't know what to use, use `search_ng21`.

| search algorithm (inside `fmindex_collection::`)                | Description |
|:----------------------------------------------------------------|-------------|
| `search_no_errors::search(index_t, query_t, cb_t)`              | searches for a perfect match |
| `search_one_error::search(index_t, query_t, cb_t)`              | search for a match with exactly one error, applying hamming distance|
| `search_pseudo::search(index_t, query_t, scheme_t, cb_t)`       | finds all alignments. uses a search scheme for efficient searching |
| `search_backtracking::search(index_t, queries_t, size_t, cb_t)` | Uses no search scheme, but naive backtracking |
| `search_backtracking_with_buffers::search`                      | Same as backtracking algorithm, but uses internal two buffers, is usually faster |
| `search_ng12::search(index_t, query_t, scheme_t, cb_t)`         | optimized by removing certain insert/substitution/deletion combinations |
| `search_ng14::search(index_t, query_t, scheme_t, cb_t)`         | same as search_ng12 but with small optimizations |
| `search_ng15::search(index_t, query_t, scheme_t, cb_t)`         | same as search_ng14 but search direction is predetermined (small optimization) |
| `search_ng16::search(index_t, query_t, scheme_t, cb_t)`         | combines ng15 and ng20 into the fastest search with large allowed errors (similar to columba) |
| `search_ng17::search`                                           | Same as ng16, but with a banded matrix. |
| `search_ng17ea::search`                                         | Same as ng17 but with early abort. |
| `search_ng20::search(index_t, query_t, scheme_t, cb_t)`         | using an banded alignment matrix (only works with backtracking search schemes) |
| `search_ng21::search(index_t, query_t, scheme_t, cb_t)`         | similar to search_ng14 but with optimizations also leaving out certain merge combination if different search path exists |
| `search_ng21ea::search`                                         | same as ng21 but with early abort. |
| `search_ng21V2::search`                                         | same as ng21 but slight internal changes |
| `search_ng21V3::search`                                         | same as ng21 but slight internal changes |
| `search_ng21V4::search`                                         | same as ng21 but slight internal changes |
| `search_ng21V5::search`                                         | same as ng21 but slight internal changes |
| `search_ng21V6::search`                                         | same as ng21 but slight internal changes |
| `search_ng21V7::search`                                         | same as ng21 but slight internal changes |
| `search_ng22::search(index_t, query_t, scheme_t, cb_t)`         | same as search_ng21 but actually doesn't do a search, but an alignment |
