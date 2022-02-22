#include <fmindex-collection/occtable/all.h>
#include <fmindex-collection/BiFMIndex.h>
#include <fmindex-collection/BiFMIndex32.h>
#include <catch2/catch.hpp>


TEMPLATE_TEST_CASE("checking bidirectional fm index", "[BiFMIndex]",
    fmindex_collection::occtable::bitvector::OccTable<256>,
    fmindex_collection::occtable::bitvectorPrefix::OccTable<256>,
    fmindex_collection::occtable::compact::OccTable<256>,
    fmindex_collection::occtable::compactAligned::OccTable<256>,
    fmindex_collection::occtable::compact2::OccTable<256>,
    fmindex_collection::occtable::compact2Aligned::OccTable<256>,
    fmindex_collection::occtable::compactPrefix::OccTable<256>,
    fmindex_collection::occtable::wavelet::OccTable<256>,
    fmindex_collection::occtable::compactWavelet::OccTable<256>,
    fmindex_collection::occtable::naive::OccTable<256>,
    fmindex_collection::occtable::sdsl_wt_bldc::OccTable<256>
) {
    using OccTable = TestType;

    auto bwt    = std::vector<uint8_t>{'t', '\0', 'o', '\0', ' ', 'H', 'W', 'a', 'l', 'e', 'l', 'l'};
    auto bwtRev = std::vector<uint8_t>{'H', '\0', 'W', 'a', 'e', 'l', 'l', 'l', 't', 'o', ' ', '\0'};
    auto sa     = std::vector<size_t>{ 10, 11, 5, 0,  6,  1,  7,  2,  3,  8,  4,  9 };

    SECTION("full sa") {
        auto bitStack = fmindex_collection::BitStack{};
        for (size_t i{0}; i < sa.size(); ++i) {
            bitStack.push(true);
        }
        auto csa = fmindex_collection::CSA{sa, bitStack};
        auto index = fmindex_collection::BiFMIndex<OccTable>{bwt, bwtRev, std::move(csa)};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            CHECK(index.locate(i) == std::make_tuple(0, sa[i]));
        }
    }

    SECTION("sa with only every second value given - sa sampled") {
        auto bitStack = fmindex_collection::BitStack{};
        auto sa2 = std::vector<size_t>{};
        for (size_t i{0}; i < sa.size(); ++i) {
            auto add = bool{i % 2 == 0} || (sa[i] == 0);
            bitStack.push(add);
            if (add) {
                sa2.push_back(sa[i]);
            }
        }

        auto csa = fmindex_collection::CSA{sa2, bitStack};
        auto index = fmindex_collection::BiFMIndex<OccTable>{bwt, bwtRev, std::move(csa)};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            CHECK(index.locate(i) == std::make_tuple(0, sa[i]));
        }
    }

    SECTION("sa with only every second value given - sa sampled - uneven") {
        auto bitStack = fmindex_collection::BitStack{};
        auto sa2 = std::vector<size_t>{};
        for (size_t i{0}; i < sa.size(); ++i) {
            auto add = bool{i % 2 == 1};
            bitStack.push(add);
            if (add) {
                sa2.push_back(sa[i]);
            }
        }

        auto csa = fmindex_collection::CSA{sa2, bitStack};
        auto index = fmindex_collection::BiFMIndex<OccTable>{bwt, bwtRev, std::move(csa)};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            CHECK(index.locate(i) == std::make_tuple(0, sa[i]));
        }
    }


    SECTION("sa with only every second value given - text sampled") {
        auto bitStack = fmindex_collection::BitStack{};
        auto sa2 = std::vector<size_t>{};
        for (size_t i{0}; i < sa.size(); ++i) {
            auto add = bool{sa[i] % 2 == 0};
            bitStack.push(add);
            if (add) {
                sa2.push_back(sa[i]);
            }
        }

        auto csa = fmindex_collection::CSA{sa2, bitStack};
        auto index = fmindex_collection::BiFMIndex<OccTable>{bwt, bwtRev, std::move(csa)};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            CHECK(index.locate(i) == std::make_tuple(0, sa[i]));
        }

    }
}

TEMPLATE_TEST_CASE("checking bidirectional 32bit fm index", "[BiFMIndex32]",
    fmindex_collection::occtable::bitvector::OccTable<256>/*,
    fmindex_collection::occtable::bitvectorPrefix::OccTable<256>,
    fmindex_collection::occtable::compact::OccTable<256>,
    fmindex_collection::occtable::compactAligned::OccTable<256>,
    fmindex_collection::occtable::compact2::OccTable<256>,
    fmindex_collection::occtable::compact2Aligned::OccTable<256>,
    fmindex_collection::occtable::compactPrefix::OccTable<256>,
    fmindex_collection::occtable::wavelet::OccTable<256>,
    fmindex_collection::occtable::compactWavelet::OccTable<256>,
    fmindex_collection::occtable::naive::OccTable<256>,
    fmindex_collection::occtable::sdsl_wt_bldc::OccTable<256>*/
) {
    using OccTable = TestType;

    auto bwt    = std::vector<uint8_t>{'t', '\0', 'o', '\0', ' ', 'H', 'W', 'a', 'l', 'e', 'l', 'l'};
    auto bwtRev = std::vector<uint8_t>{'H', '\0', 'W', 'a', 'e', 'l', 'l', 'l', 't', 'o', ' ', '\0'};
    auto sa     = std::vector<size_t>{ 10, 11, 5, 0,  6,  1,  7,  2,  3,  8,  4,  9 };

    SECTION("full sa") {
        auto bitStack = fmindex_collection::BitStack{};
        for (size_t i{0}; i < sa.size(); ++i) {
            bitStack.push(true);
        }
        auto csa = fmindex_collection::CSA32{sa, bitStack};
        auto index = fmindex_collection::BiFMIndex32<OccTable>{bwt, bwtRev, std::move(csa)};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            CHECK(index.locate(i) == std::make_tuple(0, sa[i]));
        }
    }

    SECTION("sa with only every second value given - sa sampled") {
        auto bitStack = fmindex_collection::BitStack{};
        auto sa2 = std::vector<size_t>{};
        for (size_t i{0}; i < sa.size(); ++i) {
            auto add = bool{i % 2 == 0} || (sa[i] == 0);
            bitStack.push(add);
            if (add) {
                sa2.push_back(sa[i]);
            }
        }

        auto csa = fmindex_collection::CSA{sa2, bitStack};
        auto index = fmindex_collection::BiFMIndex<OccTable>{bwt, bwtRev, std::move(csa)};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            CHECK(index.locate(i) == std::make_tuple(0, sa[i]));
        }
    }

    SECTION("sa with only every second value given - sa sampled - uneven") {
        auto bitStack = fmindex_collection::BitStack{};
        auto sa2 = std::vector<size_t>{};
        for (size_t i{0}; i < sa.size(); ++i) {
            auto add = bool{i % 2 == 1};
            bitStack.push(add);
            if (add) {
                sa2.push_back(sa[i]);
            }
        }

        auto csa = fmindex_collection::CSA{sa2, bitStack};
        auto index = fmindex_collection::BiFMIndex<OccTable>{bwt, bwtRev, std::move(csa)};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            CHECK(index.locate(i) == std::make_tuple(0, sa[i]));
        }
    }


    SECTION("sa with only every second value given - text sampled") {
        auto bitStack = fmindex_collection::BitStack{};
        auto sa2 = std::vector<size_t>{};
        for (size_t i{0}; i < sa.size(); ++i) {
            auto add = bool{sa[i] % 2 == 0};
            bitStack.push(add);
            if (add) {
                sa2.push_back(sa[i]);
            }
        }

        auto csa = fmindex_collection::CSA{sa2, bitStack};
        auto index = fmindex_collection::BiFMIndex<OccTable>{bwt, bwtRev, std::move(csa)};

        REQUIRE(index.size() == bwt.size());
        for (size_t i{0}; i < sa.size(); ++i) {
            CHECK(index.locate(i) == std::make_tuple(0, sa[i]));
        }

    }
}
