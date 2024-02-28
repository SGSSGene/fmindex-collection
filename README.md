<!--
    SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
    SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
    SPDX-License-Identifier: CC-BY-4.0
-->

# FMIndex-Collection

This Repository implements data structures and algorithms that involve (bidirectional) fm indices and searching.
FMIndex structures can be used to search through huge data.
Example: The human genome consist of around 3 billion base pairs, which can be indexed into a data structure of the size of 6GB (including fmindex and compressed suffix array). This can be searched with thousands of queries per second.

See [documentation](https://sgssgene.github.io/fmindex-collection/) for detailed information.
