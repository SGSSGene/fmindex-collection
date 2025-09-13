// SPDX-FileCopyrightText: 2025 Simon Gene Gottlieb
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include <fmindex-collection/string/FlattenedBitvectors_L0L1.h>

extern "C" {
#include "AwFmIndex.h"
}

#include <vector>


//!HACK: making symbols from awFmIndex library available, even though the library does not want us to have them
extern "C" {
auto awFmIndexAlloc(const struct AwFmIndexConfiguration* _RESTRICT_ const config, const size_t bwtLength) -> struct AwFmIndex*;
auto awFmGetBlockIndexFromGlobalPosition(const size_t globalQueryPosition) -> size_t;
auto awFmGetBlockQueryPositionFromGlobalPosition(const size_t globalQueryPosition) -> uint_fast8_t;
auto AwFmMaskedVectorPopcount(const AwFmSimdVec256 vec, const uint8_t localQueryPosition) -> uint32_t;
auto awFmNucleotideCompressedVectorToLetterIndex(const uint8_t compressedVectorLetter) -> uint8_t;
auto awFmMakeNucleotideOccurrenceVector(const struct AwFmNucleotideBlock* _RESTRICT_ const blockPtr, const uint8_t letter) -> AwFmSimdVec256;
auto awFmNucleotideLetterIndexToCompressedVector(const uint8_t asciiLetter) -> uint8_t;
auto awFmMakeAminoAcidOccurrenceVector(const struct AwFmAminoBlock* _RESTRICT_ const blockPtr, const uint8_t letter) -> AwFmSimdVec256;
auto awFmAminoAcidLetterIndexToCompressedVector(const uint8_t asciiLetter) -> uint8_t;
auto awFmAminoAcidCompressedVectorToLetterIndex(const uint8_t compressedVectorLetter) -> uint8_t;
}

//!HACK: a similar function exists in AwFmIndex under the name "setBwtAndPrefixSums", but requires extra information like the suffix array.
//       This variant reduces it down to the string with rank support building part.
inline void awFmIndexCreateStringWithRankSupport(struct AwFmIndex* _RESTRICT_ const index, std::span<uint8_t const> _symbols) {
    if (index->config.alphabetType != AwFmAlphabetAmino) {
        constexpr size_t Sigma = 8;
        auto baseOccurrences = std::array<uint64_t, Sigma>{};
        for (uint64_t pos{}; pos < _symbols.size(); ++pos) {
            const size_t blockIndex       = pos / AW_FM_POSITIONS_PER_FM_BLOCK;
            const uint8_t positionInBlock = pos % AW_FM_POSITIONS_PER_FM_BLOCK;
            const uint8_t byteInVector    = positionInBlock / 8;
            const uint8_t bitInVectorByte = positionInBlock % 8;
            auto blockPtr                 = index->bwtBlockList.asNucleotide + blockIndex;
            auto letterBitVectorBytes     = (uint8_t *)blockPtr->letterBitVectors;

            if (__builtin_expect(positionInBlock == 0, 0)) {
                // when we start a new block, copy over the base occurrences, and
                // initialize the bit vectors while we only use 5 elements, copy over
                // all 8 (to preserve padding and so valgrind doesn't complain about
                // invalid writes)
                memcpy(blockPtr->baseOccurrences, baseOccurrences.data(), Sigma * sizeof(uint64_t));
                memset(blockPtr->letterBitVectors, 0, sizeof(AwFmSimdVec256) * AW_FM_NUCLEOTIDE_VECTORS_PER_WINDOW);
            }
            auto symb = _symbols[pos];

            uint8_t letterAsCompressedVector = awFmNucleotideLetterIndexToCompressedVector(symb);
            baseOccurrences[symb] += 1;
            letterBitVectorBytes[byteInVector +  0] |= ((letterAsCompressedVector >> 0) & 0x1) << bitInVectorByte;
            letterBitVectorBytes[byteInVector + 32] |= ((letterAsCompressedVector >> 1) & 0x1) << bitInVectorByte;
            letterBitVectorBytes[byteInVector + 64] |= ((letterAsCompressedVector >> 2) & 0x1) << bitInVectorByte;
        }

        // compute baseOccurrences
        for (uint8_t i = 1; i < Sigma; i++) {
            baseOccurrences[i] += baseOccurrences[i - 1];
        }
    } else {
        constexpr size_t Sigma = 24;
        auto baseOccurrences = std::array<uint64_t, Sigma>{};
        for (uint64_t pos{}; pos < _symbols.size(); ++pos) {
            const size_t blockIndex       = pos / AW_FM_POSITIONS_PER_FM_BLOCK;
            const uint8_t positionInBlock = pos % AW_FM_POSITIONS_PER_FM_BLOCK;
            const uint8_t byteInVector    = positionInBlock / 8;
            const uint8_t bitInVectorByte = positionInBlock % 8;
            auto blockPtr                 = index->bwtBlockList.asAmino + blockIndex;
            auto letterBitVectorBytes     = (uint8_t *)blockPtr->letterBitVectors;

            if (__builtin_expect(positionInBlock == 0, 0)) {
                // when we start a new block, copy over the base occurrences, and
                // initialize the bit vectors while we only use 5 elements, copy over
                // all 8 (to preserve padding and so valgrind doesn't complain about
                // invalid writes)
                memcpy(blockPtr->baseOccurrences, baseOccurrences.data(), Sigma * sizeof(uint64_t));
                memset(blockPtr->letterBitVectors, 0, sizeof(AwFmSimdVec256) * AW_FM_AMINO_VECTORS_PER_WINDOW);
            }
            auto symb = _symbols[pos];

            uint8_t letterAsCompressedVector = awFmAminoAcidLetterIndexToCompressedVector(symb);
            baseOccurrences[symb] += 1;
            letterBitVectorBytes[byteInVector + 0]   |= ((letterAsCompressedVector >> 0) & 0x1) << bitInVectorByte;
            letterBitVectorBytes[byteInVector + 32]  |= ((letterAsCompressedVector >> 1) & 0x1) << bitInVectorByte;
            letterBitVectorBytes[byteInVector + 64]  |= ((letterAsCompressedVector >> 2) & 0x1) << bitInVectorByte;
            letterBitVectorBytes[byteInVector + 96]  |= ((letterAsCompressedVector >> 3) & 0x1) << bitInVectorByte;
            letterBitVectorBytes[byteInVector + 128] |= ((letterAsCompressedVector >> 4) & 0x1) << bitInVectorByte;
        }

        // compute baseOccurrences
        for (uint8_t i = 1; i < Sigma; i++) {
            baseOccurrences[i] += baseOccurrences[i - 1];
        }
    }
}

// If AWFMIndex doesn't implement the correct version, fall back to FlattenedBitvectors_L0L1 version, which has similar memory usage
template <size_t TSigma>
struct AWFMIndex : fmc::string::FlattenedBitvectors_L0L1<TSigma, 64, 65536> {
    using fmc::string::FlattenedBitvectors_L0L1<TSigma, 64, 65536>::FlattenedBitvectors_L0L1;
};
/**
 * This is a quick and dirty implementation
 */

template <>
struct AWFMIndex<5> {
    static constexpr size_t Sigma = 5;

    AwFmIndex* index{};
    size_t totalLength{};

    AWFMIndex() = default;
    ~AWFMIndex() {
        if (index != nullptr) {
            //!WORKAROUND !HACK awFmDeallocIndex closes a file handle, even if it didn't open one
            index->fileHandle = fopen("/dev/zero", "r");
            awFmDeallocIndex(index);
        }
    }
    AWFMIndex(AWFMIndex const&) = delete;
    AWFMIndex(AWFMIndex&&) = delete;

    AWFMIndex(std::span<uint8_t const> _symbols) {
        auto config = AwFmIndexConfiguration{};
        config.alphabetType = AwFmAlphabetType::AwFmAlphabetDna;
        index = awFmIndexAlloc(&config, _symbols.size());

        awFmIndexCreateStringWithRankSupport(index, _symbols);
        totalLength = _symbols.size();
    }

    auto operator=(AWFMIndex const&) -> AWFMIndex& = delete;
    auto operator=(AWFMIndex&&) -> AWFMIndex& = delete;

    size_t size() const noexcept {
        return totalLength;
    }

    uint8_t symbol(size_t idx) const noexcept {
        auto pos = idx;
        uint64_t blockIndex = pos / AW_FM_POSITIONS_PER_FM_BLOCK;
        uint8_t  bitIndex   = pos % AW_FM_POSITIONS_PER_FM_BLOCK;

        auto blockPtr = index->bwtBlockList.asNucleotide + blockIndex;
        const uint8_t byteInBlock    = bitIndex / 8;
        const uint8_t bitInBlockByte = bitIndex % 8;

        const uint8_t *_RESTRICT_ const letterBytePointer = &((uint8_t *)&blockPtr->letterBitVectors)[byteInBlock];
        const uint8_t letterAsCompressedVector =  ((letterBytePointer[0]  >> bitInBlockByte) & 1) << 0
                                                | ((letterBytePointer[32] >> bitInBlockByte) & 1) << 1
                                                | ((letterBytePointer[64] >> bitInBlockByte) & 1) << 2;
        return awFmNucleotideCompressedVectorToLetterIndex(letterAsCompressedVector);
    }

    uint64_t rank(size_t idx, uint8_t symb) const noexcept {
        if (idx == 0) return 0;
        auto pos = idx-1;
        uint64_t blockIndex = pos / AW_FM_POSITIONS_PER_FM_BLOCK;
        uint8_t  bitIndex   = pos % AW_FM_POSITIONS_PER_FM_BLOCK;

        auto blockPtr         = index->bwtBlockList.asNucleotide + blockIndex;
        auto baseOccurrence   = blockPtr->baseOccurrences[symb];
        auto occurrenceVector = awFmMakeNucleotideOccurrenceVector(blockPtr, symb);
        auto vectorPopcount   = AwFmMaskedVectorPopcount(occurrenceVector, bitIndex);

        return vectorPopcount + baseOccurrence;
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
        //ar(occ, totalLength)
    }

    auto space_usage() const -> size_t {
        return (sizeof(AwFmNucleotideBlock) * totalLength) / AW_FM_POSITIONS_PER_FM_BLOCK;
    }
};

#define AWFMINDEX_TEMPLATE(BASESIGMA, SIGMA) \
    template <> struct AWFMIndex<SIGMA> : AWFMIndex<BASESIGMA> { \
        static constexpr size_t Sigma = SIGMA; \
        using AWFMIndex<BASESIGMA>::AWFMIndex; \
    };


AWFMINDEX_TEMPLATE(5, 2);
AWFMINDEX_TEMPLATE(5, 3);
AWFMINDEX_TEMPLATE(5, 4);

template <>
struct AWFMIndex<21> {
    static constexpr size_t Sigma = 21;

    AwFmIndex* index{};
    size_t totalLength{};

    AWFMIndex() = default;
    ~AWFMIndex() {
        if (index != nullptr) {
            //!WORKAROUND !HACK awFmDeallocIndex closes a file handle, even if it didn't open one
            index->fileHandle = fopen("/dev/zero", "r");
            awFmDeallocIndex(index);
        }
    }
    AWFMIndex(AWFMIndex const&) = delete;
    AWFMIndex(AWFMIndex&&) = delete;

    AWFMIndex(std::span<uint8_t const> _symbols) {
        auto config = AwFmIndexConfiguration{};
        config.alphabetType = AwFmAlphabetType::AwFmAlphabetAmino;
        index = awFmIndexAlloc(&config, _symbols.size());

        awFmIndexCreateStringWithRankSupport(index, _symbols);
        totalLength = _symbols.size();
    }

    auto operator=(AWFMIndex const&) -> AWFMIndex& = delete;
    auto operator=(AWFMIndex&&) -> AWFMIndex& = delete;

    size_t size() const noexcept {
        return totalLength;
    }

    uint8_t symbol(size_t idx) const noexcept {
        auto pos = idx;
        uint64_t blockIndex = pos / AW_FM_POSITIONS_PER_FM_BLOCK;
        uint8_t  bitIndex   = pos % AW_FM_POSITIONS_PER_FM_BLOCK;

        auto blockPtr = index->bwtBlockList.asAmino + blockIndex;
        const uint8_t byteInBlock    = bitIndex / 8;
        const uint8_t bitInBlockByte = bitIndex % 8;

        const uint8_t *_RESTRICT_ const letterBytePointer = &((uint8_t *)&blockPtr->letterBitVectors)[byteInBlock];
        const uint8_t letterAsCompressedVector =  ((letterBytePointer[0]   >> bitInBlockByte) & 1) << 0
                                                | ((letterBytePointer[32]  >> bitInBlockByte) & 1) << 1
                                                | ((letterBytePointer[64]  >> bitInBlockByte) & 1) << 2
                                                | ((letterBytePointer[96]  >> bitInBlockByte) & 1) << 3
                                                | ((letterBytePointer[128] >> bitInBlockByte) & 1) << 4;
        return awFmAminoAcidCompressedVectorToLetterIndex(letterAsCompressedVector);
    }

    uint64_t rank(size_t idx, uint8_t symb) const noexcept {
        if (idx == 0) return 0;
        auto pos = idx-1;
        uint64_t blockIndex = pos / AW_FM_POSITIONS_PER_FM_BLOCK;
        uint8_t  bitIndex   = pos % AW_FM_POSITIONS_PER_FM_BLOCK;

        auto blockPtr         = index->bwtBlockList.asAmino + blockIndex;
        auto baseOccurrence   = blockPtr->baseOccurrences[symb];
        auto occurrenceVector = awFmMakeAminoAcidOccurrenceVector(blockPtr, symb);
        auto vectorPopcount   = AwFmMaskedVectorPopcount(occurrenceVector, bitIndex);

        return vectorPopcount + baseOccurrence;
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
        //ar(occ, totalLength)
    }

    auto space_usage() const -> size_t {
        return (sizeof(AwFmAminoBlock) * totalLength) / AW_FM_POSITIONS_PER_FM_BLOCK;
    }
};
AWFMINDEX_TEMPLATE(21,  6);
AWFMINDEX_TEMPLATE(21,  7);
AWFMINDEX_TEMPLATE(21,  8);
AWFMINDEX_TEMPLATE(21,  9);
AWFMINDEX_TEMPLATE(21, 10);
AWFMINDEX_TEMPLATE(21, 11);
AWFMINDEX_TEMPLATE(21, 12);
AWFMINDEX_TEMPLATE(21, 13);
AWFMINDEX_TEMPLATE(21, 14);
AWFMINDEX_TEMPLATE(21, 15);
AWFMINDEX_TEMPLATE(21, 16);
AWFMINDEX_TEMPLATE(21, 17);
AWFMINDEX_TEMPLATE(21, 18);
AWFMINDEX_TEMPLATE(21, 19);
AWFMINDEX_TEMPLATE(21, 20);

//static_assert(checkString_c<AWFMIndex>);
