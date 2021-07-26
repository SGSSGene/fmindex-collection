#pragma once

#include "concepts.h"

#include <array>
#include <bitset>
#include <cstdint>
#include <vector>

#include <iostream>

namespace compactwavelet {

constexpr inline size_t bits_count(size_t y) {
    if (y == 0) return 1;
    size_t i{0};
    while (y != 0) {
        y = y >> 1;
        ++i;
    }
    return i;
}
constexpr inline size_t pow(size_t b, size_t y) {
    if (y == 0) return 1;
    return pow(b, (y-1)) * b;
}


template <size_t TSigma>
struct Bitvector {
    static constexpr auto bitct = bits_count(TSigma);
    static constexpr auto bvct = pow(2, bitct);

    struct Block {

        std::array<uint32_t, TSigma> blocks{};
        std::array<uint64_t, bitct>  bits{};

        uint64_t rank(uint8_t symb, uint8_t idx) const {
            auto which_bv = [](uint8_t symb, auto cb) {
                size_t id{0};
                uint8_t mask = 1u << (bitct-1);
                for (size_t b{0}; b < bitct; ++b) {
                    auto bit = (symb & mask) != 0;
                    cb(id, bit);
                    id = id * 2 +1 + bit;
                    mask = mask >> 1;
                }
            };

            size_t s = 0;
            size_t l1 = 64;
            size_t l2 = idx+1;
//            std::cout << "symb: " << int(symb) << " " << int(idx) << "\n";
            auto const& bbits = bits;
            size_t depth = 0;
            which_bv(symb, [&](size_t id, size_t bit) {
                auto bits = std::bitset<64>(bbits[depth]) >> s;
                auto c1 = (bits << (64-l1)).count();
                auto d1 = (bits << (64-l2)).count();


                if (bit == 0) {
                    s = s;
                    l2 = l2 - d1;
                    l1 = l1 - c1;
                } else {
                    s = s + l1 - c1;
                    l2 = d1;
                    l1 = c1;
                }
                depth += 1;
            });
//            std::cout << "final: " << l2 << "\n";
            return l2 + blocks[symb];
        }

        uint64_t prefix_rank(uint8_t symb, uint8_t idx) const {
            auto which_bv = [](uint8_t symb, auto cb) {
                size_t id{0};
                uint8_t mask = 1u << (bitct-1);
                for (size_t b{0}; b < bitct; ++b) {
                    auto bit = (symb & mask) != 0;
                    cb(id, bit);
                    id = id * 2 +1 + bit;
                    mask = mask >> 1;
                }
            };

            size_t s = 0;
            size_t l1 = 64;
            size_t l2 = idx+1;
            size_t a{};
            size_t depth = 0;
            auto const& bbits = bits;
            which_bv(symb+1, [&](size_t id, size_t bit) {
                auto bits = std::bitset<64>(bbits[depth]) >> s;
                auto c1 = (bits << (64-l1)).count();
                auto d1 = (bits << (64-l2)).count();

                auto c0 = l1 - c1;
                auto d0 = l2 - d1;
                if (bit == 0) {
                    s = s;
                    l2 = d0;
                    l1 = c0;
                } else {
                    a += d0;
                    s = s + c0;
                    l2 = d1;
                    l1 = c1;
                }
                depth += 1;
            });
            for (size_t i{0}; i <= symb; ++i) {
                a += blocks[i];
            }
            return a;
        }
    };

    std::vector<Block> blocks;
    std::vector<std::array<uint64_t, TSigma>> superBlocks;
    std::array<uint64_t, TSigma+1> C;

    template <typename CB>
    Bitvector(size_t length, CB cb) {
        blocks.reserve(length/64+2);

        blocks.emplace_back();
        superBlocks.emplace_back();

        std::array<uint64_t, bvct> sblock_acc{0};
        std::array<uint32_t, bvct> block_acc{0};
        std::array<uint64_t, bvct> bv_size{};
        for (auto& c : bv_size) {
            c = 1;
        }


        auto which_bv = [](uint8_t symb, auto cb) {
            size_t id{0};
            uint8_t mask = 1u << (bitct-1);
            for (size_t b{0}; b < bitct; ++b) {
                auto bit = (symb & mask) != 0;
                cb(id, bit);
                id = id * 2 +1 + bit;
                mask = mask >> 1;
            }
        };

        std::array<std::vector<uint8_t>, bvct> waveletcount{};
        auto insertCount = [&](size_t blockId) {
            auto& bits = blocks[blockId].bits;
            size_t idx{0};
            for (size_t i{0}; i < waveletcount.size(); ++i) {
//                std::cout << "id3: " << i << " " << waveletcount.size() << "\n";
                auto t = waveletcount.at(i);
                for (uint64_t bit : t) {
                    bits[idx / 64] = bits[idx / 64] | (bit << (idx % 64));
                    ++idx;
                }
            }
        };
        which_bv(bvct-1, [&](size_t id, size_t bit) {
            waveletcount[id].push_back(bit);
        });
        sblock_acc[bvct-1] += 1;
        block_acc[bvct-1]  += 1;

        for (size_t size{1}; size <= length; ++size) {
            if (size % (1ul<<32) == 0) { // new super block + new block
                superBlocks.emplace_back();
                for (size_t j{0}; j < TSigma; ++j) {
                    superBlocks.back()[j] = sblock_acc[j];
                }
                block_acc = {};
            }
            if (size % 64 == 0) { // new block
                insertCount(blocks.size()-1);

                blocks.emplace_back();
                for (size_t i{0}; i < TSigma; ++i) {
                    blocks.back().blocks[i] = block_acc[i];
                }
                for (auto& c : waveletcount) {
                    c.clear();
                }
            }

            auto blockId      = size >>  6;
            auto superBlockId = size >> 32;
//            auto bitId        = size &  63;

            auto symb = cb(size-1);

            which_bv(symb, [&](size_t id, size_t bit) {
//               std::cout << "id: " << id << " " << waveletcount.size() << "\n";
//               std::cout << waveletcount[id].size() << "\n";
               uint8_t bit2 = bit;
                waveletcount[id].push_back(bit2);
            });
//            std::cout << "done\n";


//            auto& bits = blocks[blockId].bits[symb];
//            bits = bits | (1ul << bitId);
            block_acc[symb] += 1;
            sblock_acc[symb] += 1;
        }
        while (waveletcount[0].size() < 64) {
            which_bv(bvct-1, [&](size_t id, size_t bit) {
//                     std::cout << "id2: " << id << "\n";
                waveletcount[id].push_back(bit);
            });
            block_acc[bvct-1] += 1;
            sblock_acc[bvct-1] += 1;
        }
        insertCount(blocks.size()-1);

        C[0] = 0;
        for (size_t i{0}; i < TSigma; ++i) {
            C[i+1] = sblock_acc[i] + C[i];
        }
    };


    uint64_t rank(uint8_t symb, uint64_t idx) const {
        auto blockId      = idx >>  6;
        auto superBlockId = idx >> 32;
        auto bitId        = idx &  63;
        return blocks[blockId].rank(symb, bitId) + superBlocks[superBlockId][symb] + C[symb];
    }

    uint64_t prefix_rank(uint8_t symb, uint64_t idx) const {
        auto blockId      = idx >>  6;
        auto superBlockId = idx >> 32;
        auto bitId        = idx &  63;
        uint64_t a={};
        for (size_t i{0}; i<= symb; ++i) {
            a += superBlocks[superBlockId][i];
        }
        return blocks[blockId].prefix_rank(symb, bitId) + a;
    }

};


template <size_t TSigma>
struct FMIndex {
    static constexpr size_t Sigma = TSigma;

    Bitvector<Sigma> bitvector;
    size_t size_{};

    static size_t expectedMemoryUsage(size_t length) {
        using Block = typename Bitvector<TSigma>::Block;
        auto blockSize = std::max(alignof(Block), sizeof(Block));

        size_t C           = sizeof(uint64_t) * (Sigma+1);
        size_t blocks      = blockSize        * (length+1) / 64;
        size_t superblocks = sizeof(uint64_t) * (length+1) / (1ul << 32);
        return C + blocks + superblocks;
    }


    FMIndex(std::vector<uint8_t> const& _bwt)
        : bitvector(_bwt.size(), [&](size_t i) -> uint8_t {
            return _bwt[i];
        })
        , size_{_bwt.size()}
    {}

    uint64_t size() const {
        return size_;
    }


    uint64_t rank(uint8_t symb, uint64_t idx) const {
        return bitvector.rank(symb, idx);
    }

    uint64_t prefix_rank(uint8_t symb, uint64_t idx) const {
        return bitvector.prefix_rank(symb, idx);
    }
};
static_assert(checkFMIndex<FMIndex>);

}
