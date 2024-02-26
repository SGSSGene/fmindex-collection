<!--
    SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
    SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
    SPDX-License-Identifier: CC-BY-4.0
-->
# Generator
These generators will generate search scheme based on principles or methods. These function are being
found under the namespace of `search_schemes::generator`.
To make them usable for the search algorithm also the function `expand` has to be called on them first.

## Implementations

| locate algorithm (inside `search_schemes::`) | Description |
| -------------------------------------------- | ----------- |
| `generator::backtracking(N, minK, K)`        | represents the standard backtracking with errors algorithm |
| `generator::bestKnown(N, minK, K)`           | mhm, I don't remember where I got these from |
| `generator::h2(N, minK, K)`                  | a custom heuristic to create a search scheme |
| `generator::kianfar(K)`                      | lists the schemes published by kianfar |
| `generator::kucherov(N, K)`                  | lists the schemes published by kucherov |
| `generator::optimum(minK, K)`                | the optimum search schemes, as they are provided in the SeqAn3 library |
| `generator::pigeon_trivial(minK, K)`         | generates schemes based on the pigeon hole principle |
| `generator::pigeon_opt(minK, K)`             | same as `pigeon_trivial` but with optimizations by merging certain searches |
| `generator::suffix_filter(N, minK, K)`       | generates based on the suffix filter algorithm |
| `generator::zeroOnesZero_trivial(minK, K)`   | generates based on the 01\*0 lossless seeds paper |
| `generator::zeroOnesZero_opt(minK, K)`       | same as above, but merging certain searches |
