# FMIndex-Collection
Data structures to implement, improve and compare bidirectional fm indices.
FMIndex structures can be used to search through huge data.
Example: The human genome consist of around 3 billion base pairs, which can be indexed into a data structure of the size of 6GB (including fmindex and compressed suffix array). This can be searched with thousands of queries per second.

# Data structures
Searching with a bidirectional index consist out of different parts:
1. The Occurrence Table (Data structure of the bidirectional index)
2. Compressed Suffix Array (Data structure to map the entries of the index into the actual position of the original text)
3. search algorithm - there are many ways to search, edit distance, hamming distance or using search schemes
4. Search scheme generator - The currently best performing search algorithms are using search schemes, these need to be generated.

## OccTables (Occurrence Tables)
Currently following structures are available and they all fulfill the "OccTable" concept.
- Naive - storing the occ table in std::vector<size_t> tables, needs O(|Σ|·n·sizeof(size_t)) space. (144GB for the human genome)
- Bitvector - using bitvector for each table O(|Σ| · n · 2/8). (4.5GB for human genome)
- BitvectorPrefix - like bitvector, but using to compute the internal prefix ranks (9GB for human genome)
- Wavelet - Using the wavelet tree structure, needs O(log|Σ| · n · 2/8) (2GB for human genome)
- Sdsl_wt_bldc - Using the wavelet tree structure implemented in the sdsl library
- Compact - using Interleaved Bitvectors, using uint32_t as block size
- CompactAligned - same as Compact, but the blocks are also memory aligned
- Compact2 - using Interleaved Bitvectors, using uint16_t as block size
- Compact2Aligned - same as Compact2, but the blocks are also memory aligned

## Compressed Suffix Array
Currently only one implementation exists
- CSA - the compressed suffix array is compressed on an suffix array based sampling.


## Search Algorithms
- SearchNoErrors - takes a list of queries and returns the result cursors
- search_pseudo - finds all solutions that correspond to a different alignment.
- search_ng12 - optimized by removing certain insert/substitution/deletion combinations
- search_ng14 - same as search_ng12 but with small optimizations
- search_ng15 - same as search_ng14 but search direction is predetermined (small optimization)
- search_ng16 - combines ng15 and ng20 into the fastest search with large allowed errors
- search_ng20 - using an banded alignment matrix (only works with backtracking search schemes)
- search_ng21 - similar to search_ng14 but with optimizations also leaving out certain merge combination if different search path exists
- search_ng22 - same as search_ng22 but actually doesn't do a search, but an alignment

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
