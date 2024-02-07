# Occurrence Table
## OccTable (`fmindex_collection::OccTable`)
Classes fulfilling this concept provide the member functions:

- `#!c++ auto rank(size_t idx, uint8_t symb) const -> size_t`
- `#!c++ auto prefix_rank(size_t idx, uint8_t symb) const -> size_t`
- `#!c++ auto all_ranks(size_t) const -> std::tuple<std::array<size_t, T::Sigma>, std::array<size_t, T::Sigma>>`
- `#!c++ auto size() const -> size_t`

## Implementation

- `#!c++ fmindex_collection::GenericOccTable<SymbolVector Vector>`
- `#!c++ fmindex_collection::occtable::compactBitvectorPrefix::OccTable<size_t TSigma>`
- `#!c++ fmindex_collection::occtable::eprV8::OccTable<size_t TSigma>`
- `#!c++ fmindex_collection::occtable::interleavedEPRV7_32::OccTablee<size_t TSigma>`
- `#!c++ fmindex_collection::occtable::interleavedPrefix::OccTablee<size_t TSigma>`
- `#!c++ fmindex_collection::occtable::interleavedWavelet32::OccTablee<size_t TSigma>`
- `#!c++ fmindex_collection::occtable::interleavedWavelet32Aligned::OccTablee<size_t TSigma>`


