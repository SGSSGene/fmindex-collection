#pragma once

#include <array>
#include <bitset>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <vector>

struct Bitvector {
    struct Superblock {
        uint64_t superBlockEntry{};
        uint64_t blockEntries{};
        std::array<uint64_t, 6> bits{};

        uint64_t rank(size_t idx) const noexcept {
            assert(idx < 384);

            auto blockId = idx >> 6;
            auto block = 0b111111111ul & (blockEntries >> (blockId * 9));
            auto keep = (idx & 63);
            auto maskedBits = bits[blockId] << (63-keep);
            auto ct = std::bitset<64>{maskedBits}.count();

            auto total = superBlockEntry + block + ct;
            return total;
        }

        bool value(size_t idx) const noexcept {
            auto blockId = idx >> 6;
            auto bitNbr = idx & 63;
            return (bits[blockId] & (1ul << bitNbr));
        }

        void setBlock(size_t blockId, size_t value) {
            blockEntries = blockEntries & ~uint64_t{0b111111111ul << blockId*9};
            blockEntries = blockEntries | uint64_t{value << blockId*9};
        }
    };


    std::vector<Superblock> superblocks{};

    size_t memoryUsage() const {
        return sizeof(superblocks) + superblocks.size() * sizeof(superblocks.back());
    }


    uint64_t rank(size_t idx) const noexcept {
        auto superblockId = idx / 384;
        auto bitId        = idx % 384;
        return superblocks[superblockId].rank(bitId);
    }

    bool value(size_t idx) const noexcept {
        idx += 1;
        auto superblockId = idx / 384;
        auto bitId        = idx % 384;
        return superblocks[superblockId].value(bitId);
    }

    template <typename CB>
    Bitvector(size_t length, CB cb) {
        superblocks.reserve(length/384+1);
        superblocks.emplace_back();
        uint64_t sblock_acc{};
        uint16_t block_acc{};

        for (size_t size{1}; size <= length; ++size) {
            if (size % 384 == 0) { // new super block + new block
                superblocks.emplace_back();
                superblocks.back().superBlockEntry = sblock_acc;
                block_acc = 0;
            } else if (size % 64 == 0) { // new block
                superblocks.back().setBlock((size % 384) / 64, block_acc);
            }

            auto blockId      = (size >>  6) % 6;
            auto superBlockId = size / 384;
            auto bitId        = size &  63;

            if (cb(size-1)) {
                auto& bits = superblocks.back().bits[blockId];
                bits = bits | (1ul << bitId);

                block_acc  += 1;
                sblock_acc += 1;
            }
        }
    }
};

