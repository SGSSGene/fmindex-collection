#pragma once

#include "concepts.h"

#include <array>
#include <bitset>
#include <cstdint>
#include <vector>


namespace fmindex_collection {
namespace occtable {
namespace interleavedEPR_impl {

// counts how many bits are needed to represent the number y
constexpr inline uint64_t bits_count(uint64_t y) {
    if (y == 0) return 1;
    uint64_t i{0};
    while (y != 0) {
        y = y >> 1;
        ++i;
    }
    return i;
}

// computes b to the power of y
constexpr inline uint64_t pow(uint64_t b, uint64_t y) {
    if (y == 0) return 1;
    return pow(b, (y-1)) * b;
}


template <uint64_t TSigma, uint64_t TAlignment, typename block_t>
struct Bitvector {

    // number of full length bitvectors needed `2^bitct ≥ TSigma`
    static constexpr auto bitct = bits_count(TSigma-1);
    // next full power of 2
    static constexpr auto bvct  = pow(2, bitct);


    // To select a char at the even/uneven position
    static constexpr uint64_t maskEven = []() {
        uint64_t entries = 64 / bitct;
        auto result = uint64_t{0};
        auto chunkMaskEven = uint64_t{(1ul << bitct)-1ul};
        for (uint64_t i{0}; i < entries; i += 2) {
            result = (result << (bitct*2)) | chunkMaskEven;
        }
        return result;
    }();

    static constexpr uint64_t bitMask = []() {
        uint64_t entries = 64 / bitct;
        auto result = uint64_t{0};
        auto mask = 1<<bitct;
        for (uint64_t i{0}; i < entries; i += 2) {
            result = (result << (bitct*2)) | mask;
        }
        return result;

    }();

    static constexpr std::array<uint64_t, TSigma> rb = []() {
        auto results = std::array<uint64_t, TSigma>{};

        for (uint64_t symb{0}; symb < TSigma; ++symb) {
            uint64_t mask = symb | (1<<bitct);
            uint64_t entries = 64 / bitct;
            for (uint64_t i{0}; i < entries; i += 2) {
                results[symb] = (results[symb] << (bitct*2)) | mask;
            }
        }

        return results;
    }();

    struct alignas(TAlignment) Block {
        std::array<block_t, TSigma> blocks{};
        uint64_t inBlock;

        void prefetch() const {
            __builtin_prefetch(reinterpret_cast<void const*>(&blocks), 0, 0);
//            __builtin_prefetch((const void*)(&bits), 0, 0);
        }


        uint64_t prefix_rank(uint64_t idx, uint64_t symb) const {
            assert(idx < 64 / bitct);

            auto _inblock = inBlock;// & ((1ul<<(idx*bitct)) -1);

            auto te = ((rb[symb] - (_inblock & maskEven)) & bitMask) >> bitct;
            auto to = (rb[symb] - ((_inblock>>bitct) & maskEven)) & bitMask;
            auto epr = (te | to) & ((1ul << (idx*bitct))-1ul);

            auto ct = std::bitset<64>{epr}.count();

            for (uint64_t i{0}; i <= symb; ++i) {
                ct += blocks[i];
            }
            return ct;
        }


        uint64_t symbol(uint64_t idx) const {
            assert(idx < 64 / bitct);

            auto mask = uint64_t{(1ul<<bitct)-1ul};
            uint64_t symb = (inBlock >> (idx*bitct)) & mask;
            return symb;
        }

        template <typename Archive>
        void serialize(Archive& ar) {
            ar(blocks, inBlock);
        }
    };

    constexpr static uint64_t letterFit = 64 / bitct;
    static constexpr uint64_t block_size = ((1ul<<(sizeof(block_t)*8)) / letterFit)*letterFit;

    std::vector<Block> blocks;
    std::vector<std::array<uint64_t, TSigma>> superBlocks;
    std::array<uint64_t, TSigma+1> C;


    template <typename CB>
    Bitvector(uint64_t length, CB cb) {
        blocks.reserve(length/64+2);

        std::array<uint64_t, TSigma> sblock_acc{0};
        std::array<block_t, TSigma> block_acc{0};

        for (uint64_t size{0}; size < length; ++size) {
            if (size % block_size == 0) { // new super block + new block
                superBlocks.emplace_back(sblock_acc);
                blocks.emplace_back();
                block_acc = {};
            } else if (size % letterFit == 0) { // new block
                blocks.emplace_back();
                blocks.back().blocks = block_acc;
            }
            auto blockId      = size / letterFit;
            auto bitId        = size % letterFit;

            uint64_t symb = cb(size);
            blocks[blockId].inBlock |= symb << (bitct * bitId);

            block_acc[symb] += 1;
            sblock_acc[symb] += 1;
        }

        C[0] = 0;
        for (uint64_t i{0}; i < TSigma; ++i) {
            C[i+1] = sblock_acc[i] + C[i];
        }
    }

    Bitvector(cereal_tag) {}


    uint64_t memoryUsage() const {
        return blocks.size() * sizeof(blocks.back())
            + superBlocks.size() * sizeof(superBlocks.back())
            + sizeof(C);
    }

    void prefetch(uint64_t idx) const {
        auto blockId      = idx / letterFit;
        blocks[blockId].prefetch();
    }

    uint64_t rank(uint64_t idx, uint64_t symb) const {
        auto blockId      = idx / letterFit;
        auto superBlockId = idx / block_size;
        auto bitId        = idx % letterFit;
        return blocks[blockId].prefix_rank(bitId, symb)
               - ((symb>0)?(blocks[blockId].prefix_rank(bitId, symb-1)):0)
               + superBlocks[superBlockId][symb]
               + C[symb];
    }

    uint64_t prefix_rank(uint64_t idx, uint64_t symb) const {
        auto blockId      = idx / letterFit;
        auto superBlockId = idx / block_size;
        auto bitId        = idx % letterFit;
        uint64_t a={};
        for (uint64_t i{0}; i<= symb; ++i) {
            a += superBlocks[superBlockId][i];
        }
        return blocks[blockId].prefix_rank(bitId, symb) + a;
    }


    auto all_ranks(uint64_t idx) const -> std::array<uint64_t, TSigma> {
        auto blockId      = idx / letterFit;
        auto superBlockId = idx / block_size;
        auto bitId        = idx % letterFit;
        auto res = std::array<uint64_t, TSigma>{};

        auto pre = std::array<uint64_t, TSigma>{};
        for (uint64_t symb{0}; symb < TSigma; ++symb) {
            pre[symb] = blocks[blockId].prefix_rank(bitId, symb);
        }

        res[0] = pre[0] + superBlocks[superBlockId][0] + C[0];

        for (uint64_t symb{1}; symb < TSigma; ++symb) {
            res[symb] = pre[symb] - pre[symb-1] + superBlocks[superBlockId][symb] + C[symb];
        }
        return res;
    }

    auto all_ranks_and_prefix_ranks(uint64_t idx) const -> std::tuple<std::array<uint64_t, TSigma>, std::array<uint64_t, TSigma>> {
        auto blockId      = idx / letterFit;
        auto superBlockId = idx / block_size;
        auto bitId        = idx % letterFit;

        auto rs  = std::array<uint64_t, TSigma>{};
        auto prs = std::array<uint64_t, TSigma>{};

        auto pre = std::array<uint64_t, TSigma>{};
        for (uint64_t symb{0}; symb < TSigma; ++symb) {
            pre[symb] = blocks[blockId].prefix_rank(bitId, symb);
        }

        rs[0] = pre[0];
        for (uint64_t symb{1}; symb < TSigma; ++symb) {
            rs[symb]  = pre[symb]-pre[symb-1];
        }
        for (uint64_t symb{0}; symb < TSigma; ++symb) {
            rs[symb] += superBlocks[superBlockId][symb];
        }
        prs[0] = rs[0];
        for (uint64_t symb{1}; symb < TSigma; ++symb) {
            prs[symb]= prs[symb-1] + rs[symb];
        }

        for (uint64_t symb{0}; symb < TSigma; ++symb) {
            rs[symb] += C[symb];
        }
        return {rs, prs};
    }

    uint64_t symbol(uint64_t idx) const {
        auto blockId      = idx / letterFit;
        auto bitId        = idx % letterFit;
        return blocks[blockId].symbol(bitId);
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(blocks, superBlocks, C);
    }
};


template <uint64_t TSigma, typename block_t, uint64_t TAlignment>
struct OccTable {
    using TLengthType = uint64_t;
    static constexpr uint64_t Sigma = TSigma;

    Bitvector<Sigma, TAlignment, block_t> bitvector;

    static uint64_t expectedMemoryUsage(uint64_t length) {
        using Block = typename Bitvector<TSigma, TAlignment, block_t>::Block;
        auto blockSize = std::max(alignof(Block), sizeof(Block));

        uint64_t C           = sizeof(uint64_t) * (Sigma+1);
        uint64_t blocks      = blockSize        * (length+1) / 64;
        uint64_t superblocks = sizeof(uint64_t) * (length+1) / (1ul << (sizeof(block_t) * 8));
        return C + blocks + superblocks;
    }

    OccTable(std::vector<uint8_t> const& _bwt)
        : bitvector(_bwt.size(), [&](uint64_t i) -> uint8_t {
            return _bwt[i];
        })
    {}

    OccTable(cereal_tag)
        : bitvector(cereal_tag{})
    {}

    uint64_t memoryUsage() const {
        return bitvector.memoryUsage() + sizeof(OccTable);
    }

    uint64_t size() const {
        return bitvector.C.back();
    }

    auto prefetch(uint64_t idx) const {
        bitvector.prefetch(idx);
    }

    uint64_t rank(uint64_t idx, uint64_t symb) const {
        return bitvector.rank(idx, symb);
    }

    uint64_t prefix_rank(uint64_t idx, uint64_t symb) const {
        return bitvector.prefix_rank(idx, symb);
    }

    uint64_t symbol(uint64_t idx) const {
        return bitvector.symbol(idx);
    }

    auto all_ranks(uint64_t idx) const -> std::tuple<std::array<uint64_t, Sigma>, std::array<uint64_t, Sigma>> {
        auto [rs, prs] = bitvector.all_ranks_and_prefix_ranks(idx);
        return {rs, prs};
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(bitvector);
    }
};
}

namespace interleavedEPR8 {
template <uint64_t TSigma>
struct OccTable : interleavedEPR_impl::OccTable<TSigma, uint8_t, 1> {
    static auto name() -> std::string {
        return "Interleaved EPR (8bit)";
    }

    static auto extension() -> std::string {
        return "iepr8";
    }
};
static_assert(checkOccTable<OccTable>);
}
namespace interleavedEPR16 {
template <uint64_t TSigma>
struct OccTable : interleavedEPR_impl::OccTable<TSigma, uint16_t, 1> {
    static auto name() -> std::string {
        return "Interleaved EPR (16bit)";
    }

    static auto extension() -> std::string {
        return "iepr16";
    }
};
static_assert(checkOccTable<OccTable>);
}
namespace interleavedEPR32 {
template <uint64_t TSigma>
struct OccTable : interleavedEPR_impl::OccTable<TSigma, uint32_t, 1> {
    static auto name() -> std::string {
        return "Interleaved EPR (32bit)";
    }

    static auto extension() -> std::string {
        return "iepr32";
    }
};
static_assert(checkOccTable<OccTable>);
}
namespace interleavedEPR8Aligned {
template <uint64_t TSigma>
struct OccTable : interleavedEPR_impl::OccTable<TSigma, uint8_t, 64> {
    static auto name() -> std::string {
        return "Interleaved EPR (8bit, aligned)";
    }

    static auto extension() -> std::string {
        return "iepr8a";
    }
};
static_assert(checkOccTable<OccTable>);
}
namespace interleavedEPR16Aligned {
template <uint64_t TSigma>
struct OccTable : interleavedEPR_impl::OccTable<TSigma, uint16_t, 64> {
    static auto name() -> std::string {
        return "Interleaved EPR (16bit, aligned)";
    }

    static auto extension() -> std::string {
        return "iepr16a";
    }
};
static_assert(checkOccTable<OccTable>);
}
namespace interleavedEPR32Aligned {
template <uint64_t TSigma>
struct OccTable : interleavedEPR_impl::OccTable<TSigma, uint32_t, 64> {
    static auto name() -> std::string {
        return "Interleaved EPR (32bit, aligned)";
    }

    static auto extension() -> std::string {
        return "iepr32a";
    }
};
static_assert(checkOccTable<OccTable>);

}

}
}
