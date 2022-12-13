# FMIndex-Collection

This Repository implements data structures and algorithms that involve (bidirectional) fm indices and searching.
FMIndex structures can be used to search through huge data.
Example: The human genome consist of around 3 billion base pairs, which can be indexed into a data structure of the size of 6GB (including fmindex and compressed suffix array). This can be searched with thousands of queries per second.

# Data structures and Algorithms
Searching with a bidirectional index consist out of different parts:
1. The **Occurrence Table** - Data structure of the bidirectional index
2. (Compressed) **Suffix Array** - Data structure to map the entries of the index into the actual position of the original text
3. **FMIndex** and **BiFMIndex** - Data structure that aggregate Occurence Table, Suffix Array and the C array into one container.
3. **Search algorithm** - there are many ways to search, edit distance, hamming distance or using search schemes
4. **Search scheme generator** - The currently best performing search algorithms are using search schemes, these need to be generated.

## OccTables (Occurrence Tables)
Currently, following structures are available and they all fulfill the "OccTable" concept. All are located in the namespace ``fmindex_collection::occtable``;
| Name                         | Type                                       | Description |
| ------------------           | -------------------------------------      | ----------- |
| **Naive**                    | `naive::OccTable<uint64_t>`                | storing the occ table in std::vector<size_t> tables, needs O(|Σ|·n·sizeof(size_t)) space. (144GB for the human genome)|
| **Bitvector**                | `bitvector::OccTable<uint64_t>`            | using bitvector for each table O(|Σ| · n · 2/8). (4.5GB for human genome) |
| **Interleaved-8**            | `interleaved8::OccTable<uint64_t>`         | using bitvectors, but interleaving the bitvectors of different occ columns |
| **Interleaved-16**           | `interleaved16::OccTable<uint64_t>`        | using bitvectors, but interleaving the bitvectors of different occ columns |
| **Interleaved-32**           | `interleaved32::OccTable<uint64_t>`        | using bitvectors, but interleaving the bitvectors of different occ columns |
| **Prefix**                   | `bitvectorPrefix::OccTable<uint64_t>`      | like bitvector, but using to compute the internal prefix ranks (9GB for human genome) |
| **Wavelet Trees**            | `wavelet::OccTable<uint64_t>`              | Using the wavelet tree structure, needs O(log|Σ| · n · 2/8) (2GB for human genome) |
| **Interleaved Wavelets**     | `interleavedWavelet::OccTable<uint64_t>`   |  |
| **Interleaved Wavelets-32**  | `interleavedWavelet32::OccTable<uint64_t>` |  |
| **Interleaved Prefix**       | `interleavedPrefix::OccTable<uint64_t>`    |  |
| **SDSL Wavelets**            | `sdsl_wt_bldc::OccTable<uint64_t>`         | Wavelets Trees based on the SDSL implementation |
| **SDSL EPR**                 | `sdsl_wt_epr::OccTable<uint64_t>`          | EPR Dictionary based on the SDSL implementation |
| **Interleaved EPR**          | `interleavedEPR::OccTable<uint64_t>`       | EPR implementation with interleaved occ tables |
| **Interleaved EPR - V2**     | `interleavedEPR::OccTable<uint64_t>`       | similar to **Interleaved EPR**, but with a different encoding on bit level |
| **EPRV3-8**                  | `epr8V3::OccTable<uint64_t>`               | similar to **Interleaved EPR V2**, but not using interleaved bitvectors, 8bit blocks |
| **EPRV3-16**                 | `epr16V3::OccTable<uint64_t>`              | similar to **Interleaved EPR V2**, but not using interleaved bitvectors, 16bit blocks |
| **EPRV3-32**                 | `epr32V3::OccTable<uint64_t>`              | similar to **Interleaved EPR V2**, but not using interleaved bitvectors, 32bit blocks |
| **EPRV4**                    | `eprV4::OccTable<uint64_t>`                | similar to **Interleaved EPR V3**, but using 8bit, 16bit, 32bit and 64bit level blocks |
| **EPRV5**                    | `eprV5::OccTable<uint64_t>`                | similar to **Interleaved EPR V3**, but using 8bit, 16bit, and 64 bit level blocks  |
| **DenseEPRV6**               | `eprV6::OccTable<uint64_t>`                | similar to **Interleaved EPR V3**, but using 8bit, 16bit, and l-bit level blocks, where l is the smallest possible block size  |

see [OccTables](doc/OccTables.png) for more details on their structure

## Compressed Suffix Array
Currently only one implementation exists
- **CSA** - the compressed suffix array is compressed on an suffix array based sampling.



## Search Algorithms
Search function are all located in the `fmindex_collection` namespace.
| Function                                                   | Description |
| -------------------------------------                      | ----------- |
| `search_no_errors::search(index_t, query_t, cb_t)`         | searches for a perfect match |
| `search_one_error::search(index_t, query_t, cb_t)`         | search for a match with exactly one error, applying hamming distance|
| `search_pseudo::search(index_t, query_t, scheme_t, cb_t)`  | finds all alignments. uses a search scheme for efficient searching |
| `search_ng12::search(index_t, query_t, scheme_t, cb_t)`    | optimized by removing certain insert/substitution/deletion combinations |
| `search_ng14::search(index_t, query_t, scheme_t, cb_t)`    | same as search_ng12 but with small optimizations |
| `search_ng15::search(index_t, query_t, scheme_t, cb_t)`    | same as search_ng14 but search direction is predetermined (small optimization) |
| `search_ng16::search(index_t, query_t, scheme_t, cb_t)`    | combines ng15 and ng20 into the fastest search with large allowed errors (similar to columba) |
| `search_ng20::search(index_t, query_t, scheme_t, cb_t)`    | using an banded alignment matrix (only works with backtracking search schemes) |
| `search_ng21::search(index_t, query_t, scheme_t, cb_t)`    | similar to search_ng14 but with optimizations also leaving out certain merge combination if different search path exists |
| `search_ng22::search(index_t, query_t, scheme_t, cb_t)`    | same as search_ng21 but actually doesn't do a search, but an alignment |

## Search Scheme generator
These generators will generate search scheme based on principles or methods
- Backtracking - represents the standard backtracking with errors algorithm
- Kianfar - lists the schemes published by kianfar
- kucherov - lists the schemes published by kucherov
- pigeon_triv - generates schemes based on the pigeon hole principle
- pigeon_opt - same as pigeon_triv but with optimizations by merging certain searches
- suffix_filter - generates based on the suffix filter algorithm
- zeroOnesZero_triv - generates based on the 01*0 lossless seeds paper
- zeroOnesZero_opt - same as above, but merging certain searches
- optimum - the optimum search schemes, as they are provided in the seqan3 library
- bestKnown - mhm, I don't remember where I got these from

