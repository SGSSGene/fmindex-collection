// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include "cereal_tag.h"

#include <array>
#include <bitset>
#include <cassert>
#if __has_include(<cereal/archives/binary.hpp>)
#include <cereal/archives/binary.hpp>
#endif
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <span>
#include <vector>

namespace fmindex_collection {

struct BitvectorCompact {
    struct Superblock {
        uint64_t superBlockEntry;
        std::array<uint8_t, 4>  blocks;
        std::array<uint64_t, 4> bits;

        uint64_t rank(size_t idx) const noexcept {
            assert(idx < 256);
            auto blockId = idx >> 6;
            auto bitId   = idx & 63;
            auto maskedBits = bits[blockId] << (63-bitId);
            auto ct = std::bitset<64>{maskedBits}.count();

            auto total = superBlockEntry + blocks[blockId] + ct;
            return total;
        }

        bool value(size_t idx) const noexcept {
            auto blockId = idx >> 6;
            auto bitId   = idx & 63;
            return (bits[blockId] >> bitId) & 1;
        }

        void setBlock(size_t blockId, size_t value) {
            blocks[blockId] = value;
        }

        template <typename Archive>
        void serialize(Archive& ar) {
            ar(superBlockEntry, blocks, bits);
        }
    };


    std::vector<Superblock> superblocks{};

    size_t memoryUsage() const {
        return sizeof(superblocks) + superblocks.size() * sizeof(superblocks.back());
    }


    uint64_t rank(size_t idx) const noexcept {
        auto superblockId = idx >> 8;
        auto bitId        = idx % 256;
        return superblocks[superblockId].rank(bitId);
    }

    bool value(size_t idx) const noexcept {
        idx += 1;
        auto superblockId = idx >> 8;
        auto bitId        = idx % 256;
        return superblocks[superblockId].value(bitId);
    }

    BitvectorCompact() = default; //!TODO should not exists

    template <typename CB>
    BitvectorCompact(size_t length, CB cb) {
        // Next three lines are a reserve call, with zero initialization
        // This is required, so padding bytes will also be zero
        superblocks.resize(length/256+1);
        memset((void*)superblocks.data(), 0, superblocks.size() * sizeof(Superblock));
        superblocks.resize(0);

        superblocks.emplace_back(Superblock{});
        uint64_t sblock_acc{};
        uint16_t block_acc{};

        for (size_t size{1}; size <= length; ++size) {
            if (size % 256 == 0) { // new super block + new block
                superblocks.emplace_back(Superblock{});
                superblocks.back().superBlockEntry = sblock_acc;
                block_acc = 0;
            } else if (size % 64 == 0) { // new block
                superblocks.back().setBlock((size % 256) / 64, block_acc);
            }

            auto blockId      = (size >>  6) % 4;
            auto bitId        = size &  63;

            if (cb(size-1)) {
                auto& bits = superblocks.back().bits[blockId];
                bits = bits | (1ull << bitId);

                block_acc  += 1;
                sblock_acc += 1;
            }
        }
    }

    BitvectorCompact(cereal_tag) {}

    template <typename Archive>
    void serialize(Archive& ar) {
        // 0 version: slow path
        // 1 version: binary path (fast)
        auto version = []() {
#if __has_include(<cereal/archives/binary.hpp>)
            if constexpr (std::same_as<Archive, cereal::BinaryOutputArchive>
                            || std::same_as<Archive, cereal::BinaryInputArchive>) {
                return 1;
            }
#endif
            return 0;
        }();
        ar(version);

        if (version == 0) {
            ar(superblocks);
        } else if (version == 1) {
#if __has_include(<cereal/archives/binary.hpp>)
            if constexpr (std::same_as<Archive, cereal::BinaryOutputArchive>
                            || std::same_as<Archive, cereal::BinaryInputArchive>) {
                auto l = superblocks.size();
                ar(l);
                superblocks.resize(l);
                ar(cereal::binary_data(superblocks.data(), l * sizeof(Superblock)));
            } else
#endif
            throw std::runtime_error("fmindex-collection - BitvectorCompact was created with binary data, but this is not available in this app");
        } else {
            throw std::runtime_error("fmindex-collection - BitvectorCompact was created with legacy format - not readable by this app");
        }
    }
};

}
