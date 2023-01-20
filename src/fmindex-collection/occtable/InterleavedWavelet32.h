#pragma once

#include "concepts.h"

#include <array>
#include <bitset>
#include <cstdint>
#include <vector>

namespace fmindex_collection {
namespace occtable {
namespace interleavedWavelet32_detail {

// counts how many bits are needed to represent the number y
constexpr inline uint64_t bits_count(uint32_t y) {
    if (y == 0) return 1;
    uint32_t i{0};
    while (y != 0) {
        y = y >> 1;
        ++i;
    }
    return i;
}

// computes b to the power of y
constexpr inline uint32_t pow(uint32_t b, uint32_t y) {
    if (y == 0) return 1;
    return pow(b, (y-1)) * b;
}

/*
 * \param TSigma size of the alphabet
 * \param TAlignment byte alignment request
 */
template <uint32_t TSigma, uint32_t TAlignment>
struct Bitvector {
    /** traversers the bits of symb
     */
    template <typename CB>
    static void traverse_symb_bits (uint32_t symb, CB cb) {
        uint32_t id{0};
        uint32_t mask = 1u << (bitct-1);
        for (uint32_t b{0}; b < bitct; ++b) {
            auto bit = (symb & mask) != 0;
            cb(id, bit);
            id = id * 2 + 1 + bit;
            mask = mask >> 1ul;
        }
    }

    // number of full length bitvectors needed `2^bitct ≥ TSigma`
    static constexpr auto bitct = bits_count(TSigma);
    // next full power of 2
    static constexpr auto bvct  = pow(2, bitct);

    struct alignas(TAlignment) Block {
        std::array<uint32_t, TSigma> blocks{};
        std::array<uint64_t, bitct>  bits{};

        auto prefetch() const {
            __builtin_prefetch(reinterpret_cast<void const*>(&blocks), 0, 0);
//            __builtin_prefetch((const void*)(&bits), 0, 0);
        }

        uint32_t rank(uint32_t idx, uint32_t symb) const {
            uint32_t s = 0;
            uint32_t l1 = 64;
            uint32_t l2 = idx;
            auto const& bbits = bits;
            uint32_t depth = 0;
            traverse_symb_bits(symb, [&]([[maybe_unused]] uint32_t id, uint32_t bit) {
                auto bits = std::bitset<64>(bbits[depth]) >> s;
                auto c1 = uint32_t((bits << (64-l1)).count());
                auto d1 = uint32_t((bits << (64-l2)).count());

                if (bit == 0) {
                    //s = s;
                    l2 = l2 - d1;
                    l1 = l1 - c1;
                } else {
                    s = s + l1 - c1;
                    l2 = d1;
                    l1 = c1;
                }
                depth += 1;
            });
            return l2 + blocks[symb];
        }

        uint32_t prefix_rank(uint32_t idx, uint32_t symb) const {
            uint32_t s = 0;
            uint32_t l1 = 64;
            uint32_t l2 = idx;
            uint32_t a{};
            uint32_t depth = 0;
            auto const& bbits = bits;
            traverse_symb_bits(symb+1, [&]([[maybe_unused]] uint32_t id, uint32_t bit) {
                auto bits = std::bitset<64>(bbits[depth]) >> s;
                auto c1 = uint32_t((bits << (64-l1)).count());
                auto d1 = uint32_t((bits << (64-l2)).count());

                auto c0 = l1 - c1;
                auto d0 = l2 - d1;
                if (bit == 0) {
                    //s = s;
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
            for (uint32_t i{0}; i <= symb; ++i) {
                a += blocks[i];
            }
            return a;
        }

        template <typename CB>
        void bvAccess(uint32_t idx, CB cb) const {
            bvAccess(0, 64, idx, cb);
        }

        template <uint32_t Depth=0, uint64_t I=0, typename CB>
        void bvAccess(uint32_t s, uint32_t l1, uint32_t idx, CB cb) const {
            auto bv = std::bitset<64>(bits[Depth]) >> s;
            auto d1 = uint32_t((bv << (64-idx)).count());
            auto d0 = idx - d1;
            if constexpr (Depth < bitct-1) {
                auto c1 = uint32_t((bv << (64-l1)).count());
                auto c0 = l1 - c1;
                bvAccess<Depth+1, (I << 1)>  (s,    c0, d0, cb);
                bvAccess<Depth+1, (I << 1)+1>(s+c0, c1, d1, cb);
            } else {
                cb(d0, d1, std::integer_sequence<uint64_t, (I << 1)>{});
            }
        }
        template <uint32_t Depth, uint64_t I, typename CB>
        void bvAccess0(CB cb) const {
            if constexpr (Depth < bitct-1) {
                bvAccess0<Depth+1, (I << 1)>(cb);
                bvAccess0<Depth+1, (I << 1)+1>(cb);
            } else {
                cb(0, 0, std::integer_sequence<uint64_t, (I << 1)>{});
            }
        }

        auto all_ranks(uint32_t idx) const -> std::array<uint32_t, TSigma> {
            std::array<uint32_t, TSigma> rs{0};
            bvAccess(idx, [&]<uint64_t I>(uint32_t d0, uint32_t d1, std::integer_sequence<uint64_t, I>) {
                if constexpr (I < TSigma) {
                    rs[I]    = d0 + blocks[I];
                }
                if constexpr (I+1 < TSigma) {
                    rs[I+1]  = d1 + blocks[I+1];
                }
            });
            return rs;
        }

        template <typename Archive>
        void serialize(Archive& ar) {
            ar(blocks, bits);
        }
    };


    std::vector<Block> blocks;
    std::array<uint32_t, TSigma+1> C;

    template <typename CB>
    Bitvector(uint32_t length, CB cb) {
        blocks.reserve(length/64+2);

        std::array<uint32_t, bvct> sblock_acc{0};
        std::array<uint32_t, bvct> block_acc{0};
        std::array<uint64_t, bvct> bv_size{};
        for (auto& c : bv_size) {
            c = 0;
        }

        std::array<std::vector<uint8_t>, bvct> waveletcount{};
        auto insertCount = [&]() {
            if (!blocks.empty()) {
                auto& bits = blocks.back().bits;
                uint32_t idx{0};
                for (uint32_t i{0}; i < waveletcount.size(); ++i) {
                    auto t = waveletcount.at(i);
                    for (uint64_t bit : t) {
                        bits[idx / 64] = bits[idx / 64] | (bit << (idx % 64));
                        ++idx;
                    }
                }
            }
            assert(blocks.size() <( 1ul<<26));

            blocks.emplace_back();
            for (uint32_t i{0}; i < TSigma; ++i) {
                blocks.back().blocks[i] = block_acc[i];
            }
            for (auto& c : waveletcount) {
                c.clear();
            }
        };

        for (uint32_t size{0}; size < length; ++size) {
            if (size % 64 == 0) { // new block
               insertCount();
            }

            auto symb = cb(size);
            traverse_symb_bits(symb, [&](uint32_t id, uint32_t bit) {
                waveletcount[id].push_back(bit);
            });
            block_acc[symb] += 1;
            sblock_acc[symb] += 1;
        }
        while (waveletcount[0].size() < 64) {
            traverse_symb_bits(bvct-1, [&](uint32_t id, uint32_t bit) {
                waveletcount[id].push_back(bit);
            });
            block_acc[bvct-1] += 1;
            sblock_acc[bvct-1] += 1;
        }
        insertCount();

        C[0] = 0;
        for (uint32_t i{0}; i < TSigma; ++i) {
            C[i+1] = sblock_acc[i] + C[i];
        }
    }

    Bitvector(cereal_tag) {}


    uint64_t memoryUsage() const {
        return blocks.size() * sizeof(blocks.back())
            + sizeof(C);
    }

    auto prefetch(uint32_t idx) const {
        auto blockId      = idx >>  6;
        blocks[blockId].prefetch();
    }

    uint32_t rank(uint32_t idx, uint32_t symb) const {
        auto blockId      = idx >>  6;
        auto bitId        = idx &  63;
        return blocks[blockId].rank(bitId, symb) + C[symb];
    }

    uint32_t prefix_rank(uint32_t idx, uint32_t symb) const {
        auto blockId      = idx >>  6;
        auto bitId        = idx &  63;
        return blocks[blockId].prefix_rank(bitId, symb);
    }

    auto all_ranks(uint32_t idx) const -> std::tuple<std::array<uint32_t, TSigma>, std::array<uint32_t, TSigma>> {
        auto blockId      = idx >>  6;
        auto bitId        = idx &  63;

        auto rs = blocks[blockId].all_ranks(bitId);

        std::array<uint32_t, TSigma> prs;

        prs[0] = rs[0];
        rs[0] += C[0];
        for (uint32_t i{1}; i < TSigma; ++i) {

            prs[i] = prs[i-1] + rs[i];
            rs[i] += C[i];
        }

        return {rs, prs};
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(blocks, C);
    }
};

/* Implements the concept `OccTable`
 *
 * \param TSigma size of the alphabet
 * \param TAlignment byte alignment request
 */
template <uint32_t TSigma, uint32_t TAlignment>
struct OccTable {
    using TLengthType = uint32_t;
    static constexpr uint32_t Sigma = TSigma;

    Bitvector<Sigma, TAlignment> bitvector;
    uint64_t size_{};

    static uint64_t expectedMemoryUsage(uint32_t length) {
        using Block = typename Bitvector<TSigma, TAlignment>::Block;
        auto blockSize = std::max(alignof(Block), sizeof(Block));

        uint64_t C           = sizeof(uint32_t) * (Sigma+1);
        uint64_t blocks      = blockSize        * (length+1) / 64;
        return C + blocks;
    }


    OccTable(std::vector<uint8_t> const& _bwt)
        : bitvector(_bwt.size(), [&](uint64_t i) -> uint8_t {
            return _bwt[i];
        })
        , size_{_bwt.size()}
    {}

    OccTable(cereal_tag)
        : bitvector{cereal_tag{}}
    {}

    uint64_t memoryUsage() const {
        return bitvector.memoryUsage() + sizeof(OccTable);
    }

    uint32_t size() const {
        return size_;
    }

    auto prefetch(uint32_t idx) const {
        bitvector.prefetch(idx);
    }

    uint32_t rank(uint32_t idx, uint64_t symb) const {
        return bitvector.rank(idx, symb);
    }

    uint32_t prefix_rank(uint32_t idx, uint64_t symb) const {
        return bitvector.prefix_rank(idx, symb);
    }

    uint32_t symbol(uint32_t idx) const {
        idx += 1;
        auto [r1, pr1] = all_ranks(idx);
        auto [r2, pr2] = all_ranks(idx-1);

        for (uint32_t i{1}; i < Sigma; ++i) {
            if (r1[i] > r2[i]) {
                return i;
            }
        }
        return 0;
    }

    auto all_ranks(uint32_t idx) const -> std::tuple<std::array<uint32_t, Sigma>, std::array<uint32_t, Sigma>> {
        return bitvector.all_ranks(idx);
    }

    template <typename Archive>
    void serialize(Archive& ar) {
        ar(bitvector, size_);
    }
};
}

namespace interleavedWavelet32 {
template <uint64_t Sigma>
struct OccTable : interleavedWavelet32_detail::OccTable<Sigma, 8> {
    static auto name() -> std::string {
        return "Interleaved Wavelet (size: 32bit)";
    }

    static auto extension() -> std::string {
        return "iw32";
    }
};
static_assert(checkOccTable<OccTable>);
}

namespace interleavedWavelet32Aligned {
template <uint64_t Sigma>
struct OccTable : interleavedWavelet32_detail::OccTable<Sigma, 64> {
    static auto name() -> std::string {
        return "Interleaved Wavelet Aligned (size: 32bit)";
    }

    static auto extension() -> std::string {
        return "iw32a";
    }
};
static_assert(checkOccTable<OccTable>);
}

}
}
