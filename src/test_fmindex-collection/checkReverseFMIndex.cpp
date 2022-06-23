#include <fmindex-collection/occtable/all.h>
#include <fmindex-collection/ReverseFMIndex.h>
#include <catch2/catch.hpp>


TEMPLATE_TEST_CASE("checking reverse unidirectional fm index", "[ReverseFMIndex]",
    fmindex_collection::occtable::bitvector::OccTable<256>,
    fmindex_collection::occtable::bitvectorPrefix::OccTable<256>,
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
    fmindex_collection::occtable::naive::OccTable<256>
//    fmindex_collection::occtable::sdsl_wt_bldc::OccTable<256>
//    fmindex_collection::occtable::sdsl_wt_epr::OccTable<256>
) {
    using OccTable = TestType;

    auto bwt = std::vector<uint8_t>{'H', '\0', 'W', 'a', 'e', 'l', 'l', 'l', 't', 'o', ' ', '\0'};
    auto sa  = std::vector<uint64_t>{0, 11, 6, 1, 7, 2, 8, 3, 9, 4, 5, 10};

    SECTION("full sa") {
        auto bitStack = fmindex_collection::BitStack{};
        for (size_t i{0}; i < sa.size(); ++i) {
            bitStack.push(true);
        }
        auto csa = fmindex_collection::CSA{sa, bitStack, 1, 63};
        auto index = fmindex_collection::ReverseFMIndex<OccTable>{bwt, std::move(csa)};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            CHECK(index.locate(i) == std::make_tuple(0, sa[i]));
        }
    }

    SECTION("sa with only every second value given - sa sampled") {
        auto bitStack = fmindex_collection::BitStack{};
        auto sa2 = std::vector<uint64_t>{};
        for (size_t i{0}; i < sa.size(); ++i) {
            auto add = bool{i % 2 == 0} || (sa[i] == sa.size()-1);
            bitStack.push(add);
            if (add) {
                sa2.push_back(sa[i]);
            }
        }

        auto csa = fmindex_collection::CSA{sa2, bitStack, 2, 63};
        auto index = fmindex_collection::ReverseFMIndex<OccTable>{bwt, std::move(csa)};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            INFO(i);
            INFO(sa[i]);
            CHECK(index.locate(i) == std::make_tuple(0, sa[i]));
        }
    }

    SECTION("sa with only every second value given - sa sampled - uneven") {
        auto bitStack = fmindex_collection::BitStack{};
        auto sa2 = std::vector<uint64_t>{};
        for (size_t i{0}; i < sa.size(); ++i) {
            auto add = bool{i % 2 == 1} || (sa[i] == sa.size()-1);
            bitStack.push(add);
            if (add) {
                sa2.push_back(sa[i]);
            }
        }

        auto csa = fmindex_collection::CSA{sa2, bitStack, 2, 63};
        auto index = fmindex_collection::ReverseFMIndex<OccTable>{bwt, std::move(csa)};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            INFO(i);
            INFO(std::get<1>(index.locate(i)));
            INFO(sa[i]);
            CHECK(index.locate(i) == std::make_tuple(0, sa[i]));
        }
    }


    SECTION("sa with only every second value given - text sampled") {
        auto bitStack = fmindex_collection::BitStack{};
        auto sa2 = std::vector<uint64_t>{};
        for (size_t i{0}; i < sa.size(); ++i) {
            auto add = bool{sa[i] % 2 == 0} || (sa[i] == 11);
            bitStack.push(add);
            if (add) {
                sa2.push_back(sa[i]);
            }
        }

        auto csa = fmindex_collection::CSA{sa2, bitStack, 2, 63};
        auto index = fmindex_collection::ReverseFMIndex<OccTable>{bwt, std::move(csa)};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            INFO(i);
            INFO(std::get<1>(index.locate(i)));
            INFO(sa[i]);
            CHECK(index.locate(i) == std::make_tuple(0, sa[i]));
        }

    }

    SECTION("compare to a directly created index") {
        auto bwt = std::vector<uint8_t>{'H', 'W', 'a', 'e', 'l', 'l', 'l', 't', 'o', ' ', '\0'};
        auto sa  = std::vector<uint64_t>{0, 6, 1, 7, 2, 8, 3, 9, 4, 5, 10};

        auto text  = std::vector<uint8_t>{'H', 'a', 'l', 'l', 'o', ' ', 'W', 'e', 'l', 't'};
        auto index = fmindex_collection::ReverseFMIndex<OccTable>{std::vector<std::vector<uint8_t>>{text}, 1};

        REQUIRE(bwt.size() == index.size());
        REQUIRE(sa.size() == index.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            INFO(i);
            INFO(sa[i]);
            INFO(std::get<1>(index.locate(i)));

            CHECK(index.locate(i) == std::make_tuple(0, sa[i]));
        }
    }

    SECTION("compare to a directly created index but with sampling") {
        auto bwt = std::vector<uint8_t>{'H', 'W', 'a', 'e', 'l', 'l', 'l', 't', 'o', ' ', '\0'};
        auto sa  = std::vector<uint64_t>{0, 6, 1, 7, 2, 8, 3, 9, 4, 5, 10};

        auto text  = std::vector<uint8_t>{'H', 'a', 'l', 'l', 'o', ' ', 'W', 'e', 'l', 't'};
        auto index = fmindex_collection::ReverseFMIndex<OccTable>{std::vector<std::vector<uint8_t>>{text}, 2};

        REQUIRE(bwt.size() == index.size());
        REQUIRE(sa.size() == index.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            INFO(i);
            INFO(sa[i]);
            INFO(std::get<1>(index.locate(i)));

            CHECK(index.locate(i) == std::make_tuple(0, sa[i]));
        }
    }

}
