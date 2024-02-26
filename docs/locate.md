<!--
    SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
    SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
    SPDX-License-Identifier: CC-BY-4.0
-->
# Locate

After searching locating is required, two algorithms are implemented

| locate algorithm (inside `fmindex_collection::`)                | Description |
|:----------------------------------------------------------------|-------------|
| `LocateLinear`                                                  | Standard linear locate |
| `LocateFMTree`                                                  | FMTree locate (not faster in this implementation) |
