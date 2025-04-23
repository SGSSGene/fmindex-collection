#pragma once

extern "C" {
#include "AwFmIndex.h"
}

#include <vector>

extern "C" {
struct AwFmIndex *
awFmIndexAlloc(const struct AwFmIndexConfiguration *_RESTRICT_ const config,
               const size_t bwtLength);

void setBwtAndPrefixSums(
    struct AwFmIndex *_RESTRICT_ const index, const size_t bwtLength,
    const uint8_t *_RESTRICT_ const sequence,
    const uint64_t *_RESTRICT_ const unsampledSuffixArray);

AwFmSimdVec256 awFmMakeNucleotideOccurrenceVector(
    const struct AwFmNucleotideBlock *_RESTRICT_ const blockPtr,
    const uint8_t letter);

size_t awFmGetBlockIndexFromGlobalPosition(const size_t globalQueryPosition);

uint_fast8_t
awFmGetBlockQueryPositionFromGlobalPosition(const size_t globalQueryPosition);

uint32_t AwFmMaskedVectorPopcount(const AwFmSimdVec256 vec,
                                  const uint8_t localQueryPosition);

}

size_t awfmindex_rank(AwFmIndex const* index, size_t pos, uint8_t symb) {
    uint64_t blockIndex = pos / AW_FM_POSITIONS_PER_FM_BLOCK;
    uint8_t  bitIndex   = pos % AW_FM_POSITIONS_PER_FM_BLOCK;

    uint64_t baseOccurrence         = index->bwtBlockList.asNucleotide[blockIndex].baseOccurrences[symb];
    AwFmSimdVec256 occurrenceVector = awFmMakeNucleotideOccurrenceVector(index->bwtBlockList.asNucleotide + blockIndex, symb);
    uint_fast16_t vectorPopcount    = AwFmMaskedVectorPopcount(occurrenceVector, bitIndex);

    return vectorPopcount + baseOccurrence;
}
#if 0
int main() {
    auto config = AwFmIndexConfiguration{};
    config.alphabetType = AwFmAlphabetType::AwFmAlphabetDna;
    auto index = awFmIndexAlloc(&config, 10);

//    auto data = std::vector<uint8_t>{1, 2, 3, 2, 1, 2, 3, 2, 1, 2};
    auto data = std::vector<uint8_t>{'A', 'C', 'G', 'C', 'A', 'T', 'G', 'C', 'A', 'C'};
    auto sa   = std::vector<uint64_t>{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    setBwtAndPrefixSums(index, 10, data.data(), sa.data());

    std::cout << "0: 0 0 0 0\n";
    for (size_t i{0}; i < 10; ++i) {
        std::cout << (i+1) << ":";
        for (size_t s{0}; s < 4; ++s) {
            std::cout << " " << rank(index, i, s);
        }
        std::cout << "\n";
    }
}
#endif

/**
 * This is a quick and dirty implementation
 */
struct AWFMIndex {
    static constexpr size_t Sigma = 5;

    AwFmIndex* index{};
    size_t totalLength{};

    AWFMIndex() = default;
    ~AWFMIndex() {
        if (index == nullptr) {
            awFmDeallocIndex(index);
        }
    }
    AWFMIndex(std::span<uint8_t const> _symbols) {
        auto config = AwFmIndexConfiguration{};
        config.alphabetType = AwFmAlphabetType::AwFmAlphabetDna;
        index = awFmIndexAlloc(&config, _symbols.size());

        // must convert symbols into something awFmIndex can understand:
        auto data = std::vector<uint8_t>{};
        data.resize(_symbols.size());
        auto sa = std::vector<uint64_t>{};
        sa.resize(_symbols.size());
        for (size_t i{0}; i < _symbols.size(); ++i) {
            if (_symbols[i] == 0) data[i] = 'A';
            if (_symbols[i] == 1) data[i] = 'C';
            if (_symbols[i] == 2) data[i] = 'G';
            if (_symbols[i] == 3) data[i] = 'T';
            if (_symbols[i] == 4) data[i] = 'x';
            sa[i] = i+1;
        }
        setBwtAndPrefixSums(index, data.size(), data.data(), sa.data());
        totalLength = _symbols.size();
    }

    size_t size() const noexcept {
        return totalLength;
    }

    uint8_t symbol(size_t idx) const noexcept {
        return idx % Sigma;
/*        idx += 1;
        for (uint64_t i{0}; i < Sigma-1; ++i) {
            if (occ[idx][i] > occ[idx-1][i]) {
                return i;
            }
        }
        return Sigma-1;*/
    }

    uint64_t rank(size_t idx, uint8_t symb) const noexcept {
        if (idx == 0) return 0;
        return awfmindex_rank(index, idx-1, symb);
    }

    uint64_t prefix_rank(size_t idx, uint8_t symb) const noexcept {
        if (idx == 0) return 0;
        uint64_t a{};
        for (uint64_t i{0}; i <= symb; ++i) {
            a += rank(idx, i);
        }
        return a;
    }

    auto all_ranks(size_t idx) const -> std::array<uint64_t, Sigma> {
        auto r = std::array<uint64_t, Sigma>{};
        for (size_t i{0}; i < Sigma; ++i) {
            r[i] = rank(idx, i);
        }
        return r;
    }

    auto all_ranks_and_prefix_ranks(size_t idx) const -> std::tuple<std::array<uint64_t, Sigma>, std::array<uint64_t, Sigma>> {
        auto rs  = all_ranks(idx);
        auto prs = rs;
        for (size_t i{1}; i < prs.size(); ++i) {
            prs[i] = prs[i] + prs[i-1];
        }
        return {rs, prs};
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(totalLength);
        //ar(occ, totalLength);
    }
};
//static_assert(checkRankVector<AWFMIndex>);
