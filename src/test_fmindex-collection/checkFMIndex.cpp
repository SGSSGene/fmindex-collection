#include <fmindex-collection/occtable/all.h>
#include <fmindex-collection/FMIndex.h>
#include <catch2/catch.hpp>


TEMPLATE_TEST_CASE("checking unidirectional fm index", "[FMIndex]",
    fmindex_collection::occtable::compactBitvector::OccTable<256>/*,
    fmindex_collection::occtable::compactBitvectorPrefix::OccTable<256>,
    fmindex_collection::occtable::interleaved8::OccTable<256>,
    fmindex_collection::occtable::interleaved16::OccTable<256>,
    fmindex_collection::occtable::interleaved32::OccTable<256>,
    fmindex_collection::occtable::interleaved8Aligned::OccTable<256>,
    fmindex_collection::occtable::interleaved16Aligned::OccTable<256>,
    fmindex_collection::occtable::interleaved32Aligned::OccTable<256>,
    fmindex_collection::occtable::interleavedPrefix::OccTable<256>,
    fmindex_collection::occtable::wavelet::OccTable<256>,
    fmindex_collection::occtable::interleavedWavelet::OccTable<256>,
    fmindex_collection::occtable::interleavedWaveletAligned::OccTable<256>,
    fmindex_collection::occtable::interleavedWavelet32::OccTable<256>,
    fmindex_collection::occtable::interleavedWavelet32Aligned::OccTable<256>,
    fmindex_collection::occtable::interleavedEPR8::OccTable<256>,
    fmindex_collection::occtable::interleavedEPR16::OccTable<256>,
    fmindex_collection::occtable::interleavedEPR32::OccTable<256>,
    fmindex_collection::occtable::interleavedEPR8Aligned::OccTable<256>,
    fmindex_collection::occtable::interleavedEPR16Aligned::OccTable<256>,
    fmindex_collection::occtable::interleavedEPR32Aligned::OccTable<256>,
    fmindex_collection::occtable::interleavedEPR8V2::OccTable<256>,
    fmindex_collection::occtable::interleavedEPR16V2::OccTable<256>,
    fmindex_collection::occtable::interleavedEPR32V2::OccTable<256>,
    fmindex_collection::occtable::interleavedEPR8V2Aligned::OccTable<256>,
    fmindex_collection::occtable::interleavedEPR16V2Aligned::OccTable<256>,
    fmindex_collection::occtable::interleavedEPR32V2Aligned::OccTable<256>,
    fmindex_collection::occtable::naive::OccTable<256>,
    fmindex_collection::occtable::sdsl_wt_bldc::OccTable<256>*/
//    fmindex_collection::occtable::sdsl_wt_epr::OccTable<256>
) {
    using OccTable = TestType;

    auto bwt    = std::vector<uint8_t>{'t', '\0', 'o', '\0', ' ', 'H', 'W', 'a', 'l', 'e', 'l', 'l'};
    auto sa     = std::vector<uint64_t>{ 10, 11, 5, 0,  6,  1,  7,  2,  3,  8,  4,  9 };

    SECTION("full sa") {
        auto bitStack = fmindex_collection::BitStack{};
        for (size_t i{0}; i < sa.size(); ++i) {
            bitStack.push(true);
        }
        auto csa = fmindex_collection::CSA{sa, bitStack, 1, 63};
        auto index = fmindex_collection::FMIndex<OccTable>{bwt, std::move(csa)};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            auto [seqId, seqPos] = index.locate(i);
            INFO(i);
            CHECK(seqId == 0);
            CHECK(seqPos == sa[i]);
        }
    }

    SECTION("sa with only every second value given - sa sampled") {
        auto bitStack = fmindex_collection::BitStack{};
        auto sa2 = std::vector<uint64_t>{};
        for (size_t i{0}; i < sa.size(); ++i) {
            auto add = bool{i % 2 == 0} || (sa[i] == 0);
            bitStack.push(add);
            if (add) {
                sa2.push_back(sa[i]);
            }
        }

        auto csa = fmindex_collection::CSA{sa2, bitStack, 2, 63};
        auto index = fmindex_collection::FMIndex<OccTable>{bwt, std::move(csa)};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            CHECK(index.locate(i) == std::make_tuple(0, sa[i]));
        }
    }

    SECTION("sa with only every second value given - sa sampled - uneven") {
        auto bitStack = fmindex_collection::BitStack{};
        auto sa2 = std::vector<uint64_t>{};
        for (size_t i{0}; i < sa.size(); ++i) {
            auto add = bool{i % 2 == 1};
            bitStack.push(add);
            if (add) {
                sa2.push_back(sa[i]);
            }
        }

        auto csa = fmindex_collection::CSA{sa2, bitStack, 2, 63};
        auto index = fmindex_collection::FMIndex<OccTable>{bwt, std::move(csa)};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            CHECK(index.locate(i) == std::make_tuple(0, sa[i]));
        }
    }


    SECTION("sa with only every second value given - text sampled") {
        auto bitStack = fmindex_collection::BitStack{};
        auto sa2 = std::vector<uint64_t>{};
        for (size_t i{0}; i < sa.size(); ++i) {
            auto add = bool{sa[i] % 2 == 0};
            bitStack.push(add);
            if (add) {
                sa2.push_back(sa[i]);
            }
        }

        auto csa = fmindex_collection::CSA{sa2, bitStack, 2, 63};
        auto index = fmindex_collection::FMIndex<OccTable>{bwt, std::move(csa)};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            CHECK(index.locate(i) == std::make_tuple(0, sa[i]));
        }

    }
}
