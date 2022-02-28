#include <fmindex-collection/occtable/all.h>
#include <catch2/catch.hpp>


TEMPLATE_TEST_CASE("check if occ table is working", "[OccTable]",
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
    fmindex_collection::occtable::interleavedEPR8V2::OccTable<256>,
    fmindex_collection::occtable::interleavedEPR16V2::OccTable<256>,
    fmindex_collection::occtable::interleavedEPR32V2::OccTable<256>,
    fmindex_collection::occtable::interleavedEPR8V2Aligned::OccTable<256>,
    fmindex_collection::occtable::interleavedEPR16V2Aligned::OccTable<256>,
    fmindex_collection::occtable::interleavedEPR32V2Aligned::OccTable<256>,
    fmindex_collection::occtable::naive::OccTable<256>,
    fmindex_collection::occtable::sdsl_wt_bldc::OccTable<256>
//    fmindex_collection::occtable::sdsl_wt_epr::OccTable<256>
) {
    using OccTable = TestType;

    auto text = std::vector<uint8_t>{'H', 'a', 'l', 'l', 'o', ' ', 'W', 'e', 'l', 't'};

    auto table = OccTable{text};

    REQUIRE(table.size() == text.size());

    SECTION("check that symbol() call works") {
        for (size_t i{0}; i < text.size(); ++i) {
            INFO(i);
            CHECK(table.symbol(i) == text.at(i));
        }
    }
    CHECK(!OccTable::name().empty());
    CHECK(!OccTable::extension().empty());
    CHECK(OccTable::Sigma == 256);
    SECTION("test complete table 'H' for rank()") {
        CHECK(table.rank( 0, ' ') == 0);
        CHECK(table.rank( 1, ' ') == 0);
        CHECK(table.rank( 2, ' ') == 0);
        CHECK(table.rank( 3, ' ') == 0);
        CHECK(table.rank( 4, ' ') == 0);
        CHECK(table.rank( 5, ' ') == 0);
        CHECK(table.rank( 6, ' ') == 1);
        CHECK(table.rank( 7, ' ') == 1);
        CHECK(table.rank( 8, ' ') == 1);
        CHECK(table.rank( 9, ' ') == 1);
        CHECK(table.rank(10, ' ') == 1);

        CHECK(table.rank( 0, 'H') == 1);
        CHECK(table.rank( 1, 'H') == 2);
        CHECK(table.rank( 2, 'H') == 2);
        CHECK(table.rank( 3, 'H') == 2);
        CHECK(table.rank( 4, 'H') == 2);
        CHECK(table.rank( 5, 'H') == 2);
        CHECK(table.rank( 6, 'H') == 2);
        CHECK(table.rank( 7, 'H') == 2);
        CHECK(table.rank( 8, 'H') == 2);
        CHECK(table.rank( 9, 'H') == 2);
        CHECK(table.rank(10, 'H') == 2);

        CHECK(table.rank( 0, 'W') == 2);
        CHECK(table.rank( 1, 'W') == 2);
        CHECK(table.rank( 2, 'W') == 2);
        CHECK(table.rank( 3, 'W') == 2);
        CHECK(table.rank( 4, 'W') == 2);
        CHECK(table.rank( 5, 'W') == 2);
        CHECK(table.rank( 6, 'W') == 2);
        CHECK(table.rank( 7, 'W') == 3);
        CHECK(table.rank( 8, 'W') == 3);
        CHECK(table.rank( 9, 'W') == 3);
        CHECK(table.rank(10, 'W') == 3);

        CHECK(table.rank( 0, 'a') == 3);
        CHECK(table.rank( 1, 'a') == 3);
        CHECK(table.rank( 2, 'a') == 4);
        CHECK(table.rank( 3, 'a') == 4);
        CHECK(table.rank( 4, 'a') == 4);
        CHECK(table.rank( 5, 'a') == 4);
        CHECK(table.rank( 6, 'a') == 4);
        CHECK(table.rank( 7, 'a') == 4);
        CHECK(table.rank( 8, 'a') == 4);
        CHECK(table.rank( 9, 'a') == 4);
        CHECK(table.rank(10, 'a') == 4);

        CHECK(table.rank( 0, 'e') == 4);
        CHECK(table.rank( 1, 'e') == 4);
        CHECK(table.rank( 2, 'e') == 4);
        CHECK(table.rank( 3, 'e') == 4);
        CHECK(table.rank( 4, 'e') == 4);
        CHECK(table.rank( 5, 'e') == 4);
        CHECK(table.rank( 6, 'e') == 4);
        CHECK(table.rank( 7, 'e') == 4);
        CHECK(table.rank( 8, 'e') == 5);
        CHECK(table.rank( 9, 'e') == 5);
        CHECK(table.rank(10, 'e') == 5);

        CHECK(table.rank( 0, 'l') == 5);
        CHECK(table.rank( 1, 'l') == 5);
        CHECK(table.rank( 2, 'l') == 5);
        CHECK(table.rank( 3, 'l') == 6);
        CHECK(table.rank( 4, 'l') == 7);
        CHECK(table.rank( 5, 'l') == 7);
        CHECK(table.rank( 6, 'l') == 7);
        CHECK(table.rank( 7, 'l') == 7);
        CHECK(table.rank( 8, 'l') == 7);
        CHECK(table.rank( 9, 'l') == 8);
        CHECK(table.rank(10, 'l') == 8);

        CHECK(table.rank( 0, 'o') == 8);
        CHECK(table.rank( 1, 'o') == 8);
        CHECK(table.rank( 2, 'o') == 8);
        CHECK(table.rank( 3, 'o') == 8);
        CHECK(table.rank( 4, 'o') == 8);
        CHECK(table.rank( 5, 'o') == 9);
        CHECK(table.rank( 6, 'o') == 9);
        CHECK(table.rank( 7, 'o') == 9);
        CHECK(table.rank( 8, 'o') == 9);
        CHECK(table.rank( 9, 'o') == 9);
        CHECK(table.rank(10, 'o') == 9);

        CHECK(table.rank( 0, 't') ==  9);
        CHECK(table.rank( 1, 't') ==  9);
        CHECK(table.rank( 2, 't') ==  9);
        CHECK(table.rank( 3, 't') ==  9);
        CHECK(table.rank( 4, 't') ==  9);
        CHECK(table.rank( 5, 't') ==  9);
        CHECK(table.rank( 6, 't') ==  9);
        CHECK(table.rank( 7, 't') ==  9);
        CHECK(table.rank( 8, 't') ==  9);
        CHECK(table.rank( 9, 't') ==  9);
        CHECK(table.rank(10, 't') == 10);
    }

    SECTION("test complete table 'H' for prefix_rank()") {
        CHECK(table.prefix_rank( 0, ' ') == 0);
        CHECK(table.prefix_rank( 1, ' ') == 0);
        CHECK(table.prefix_rank( 2, ' ') == 0);
        CHECK(table.prefix_rank( 3, ' ') == 0);
        CHECK(table.prefix_rank( 4, ' ') == 0);
        CHECK(table.prefix_rank( 5, ' ') == 0);
        CHECK(table.prefix_rank( 6, ' ') == 1);
        CHECK(table.prefix_rank( 7, ' ') == 1);
        CHECK(table.prefix_rank( 8, ' ') == 1);
        CHECK(table.prefix_rank( 9, ' ') == 1);
        CHECK(table.prefix_rank(10, ' ') == 1);

        CHECK(table.prefix_rank( 0, 'H') == 0);
        CHECK(table.prefix_rank( 1, 'H') == 1);
        CHECK(table.prefix_rank( 2, 'H') == 1);
        CHECK(table.prefix_rank( 3, 'H') == 1);
        CHECK(table.prefix_rank( 4, 'H') == 1);
        CHECK(table.prefix_rank( 5, 'H') == 1);
        CHECK(table.prefix_rank( 6, 'H') == 2);
        CHECK(table.prefix_rank( 7, 'H') == 2);
        CHECK(table.prefix_rank( 8, 'H') == 2);
        CHECK(table.prefix_rank( 9, 'H') == 2);
        CHECK(table.prefix_rank(10, 'H') == 2);

        CHECK(table.prefix_rank( 0, 'W') == 0);
        CHECK(table.prefix_rank( 1, 'W') == 1);
        CHECK(table.prefix_rank( 2, 'W') == 1);
        CHECK(table.prefix_rank( 3, 'W') == 1);
        CHECK(table.prefix_rank( 4, 'W') == 1);
        CHECK(table.prefix_rank( 5, 'W') == 1);
        CHECK(table.prefix_rank( 6, 'W') == 2);
        CHECK(table.prefix_rank( 7, 'W') == 3);
        CHECK(table.prefix_rank( 8, 'W') == 3);
        CHECK(table.prefix_rank( 9, 'W') == 3);
        CHECK(table.prefix_rank(10, 'W') == 3);

        CHECK(table.prefix_rank( 0, 'a') == 0);
        CHECK(table.prefix_rank( 1, 'a') == 1);
        CHECK(table.prefix_rank( 2, 'a') == 2);
        CHECK(table.prefix_rank( 3, 'a') == 2);
        CHECK(table.prefix_rank( 4, 'a') == 2);
        CHECK(table.prefix_rank( 5, 'a') == 2);
        CHECK(table.prefix_rank( 6, 'a') == 3);
        CHECK(table.prefix_rank( 7, 'a') == 4);
        CHECK(table.prefix_rank( 8, 'a') == 4);
        CHECK(table.prefix_rank( 9, 'a') == 4);
        CHECK(table.prefix_rank(10, 'a') == 4);

        CHECK(table.prefix_rank( 0, 'e') == 0);
        CHECK(table.prefix_rank( 1, 'e') == 1);
        CHECK(table.prefix_rank( 2, 'e') == 2);
        CHECK(table.prefix_rank( 3, 'e') == 2);
        CHECK(table.prefix_rank( 4, 'e') == 2);
        CHECK(table.prefix_rank( 5, 'e') == 2);
        CHECK(table.prefix_rank( 6, 'e') == 3);
        CHECK(table.prefix_rank( 7, 'e') == 4);
        CHECK(table.prefix_rank( 8, 'e') == 5);
        CHECK(table.prefix_rank( 9, 'e') == 5);
        CHECK(table.prefix_rank(10, 'e') == 5);

        CHECK(table.prefix_rank( 0, 'l') == 0);
        CHECK(table.prefix_rank( 1, 'l') == 1);
        CHECK(table.prefix_rank( 2, 'l') == 2);
        CHECK(table.prefix_rank( 3, 'l') == 3);
        CHECK(table.prefix_rank( 4, 'l') == 4);
        CHECK(table.prefix_rank( 5, 'l') == 4);
        CHECK(table.prefix_rank( 6, 'l') == 5);
        CHECK(table.prefix_rank( 7, 'l') == 6);
        CHECK(table.prefix_rank( 8, 'l') == 7);
        CHECK(table.prefix_rank( 9, 'l') == 8);
        CHECK(table.prefix_rank(10, 'l') == 8);

        CHECK(table.prefix_rank( 0, 'o') == 0);
        CHECK(table.prefix_rank( 1, 'o') == 1);
        CHECK(table.prefix_rank( 2, 'o') == 2);
        CHECK(table.prefix_rank( 3, 'o') == 3);
        CHECK(table.prefix_rank( 4, 'o') == 4);
        CHECK(table.prefix_rank( 5, 'o') == 5);
        CHECK(table.prefix_rank( 6, 'o') == 6);
        CHECK(table.prefix_rank( 7, 'o') == 7);
        CHECK(table.prefix_rank( 8, 'o') == 8);
        CHECK(table.prefix_rank( 9, 'o') == 9);
        CHECK(table.prefix_rank(10, 'o') == 9);

        CHECK(table.prefix_rank( 0, 't') ==  0);
        CHECK(table.prefix_rank( 1, 't') ==  1);
        CHECK(table.prefix_rank( 2, 't') ==  2);
        CHECK(table.prefix_rank( 3, 't') ==  3);
        CHECK(table.prefix_rank( 4, 't') ==  4);
        CHECK(table.prefix_rank( 5, 't') ==  5);
        CHECK(table.prefix_rank( 6, 't') ==  6);
        CHECK(table.prefix_rank( 7, 't') ==  7);
        CHECK(table.prefix_rank( 8, 't') ==  8);
        CHECK(table.prefix_rank( 9, 't') ==  9);
        CHECK(table.prefix_rank(10, 't') == 10);
    }

    SECTION("check all_ranks() is equal to prefix_rank() and rank()") {
        for (size_t idx{0}; idx < table.size(); ++idx) {
            auto [rank, prefix] = table.all_ranks(idx);
            for (size_t symb{1}; symb < 256; ++symb) {
                INFO(idx);
                INFO(symb);
                CHECK(rank[symb] == table.rank(idx, symb));
                CHECK(prefix[symb] == table.prefix_rank(idx, symb));
            }
        }
    }
}


TEMPLATE_TEST_CASE("check occ table construction on text longer than 256 characters", "[OccTable]",
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
    fmindex_collection::occtable::interleavedEPR8V2::OccTable<256>,
    fmindex_collection::occtable::interleavedEPR16V2::OccTable<256>,
    fmindex_collection::occtable::interleavedEPR32V2::OccTable<256>,
    fmindex_collection::occtable::interleavedEPR8V2Aligned::OccTable<256>,
    fmindex_collection::occtable::interleavedEPR16V2Aligned::OccTable<256>,
    fmindex_collection::occtable::interleavedEPR32V2Aligned::OccTable<256>,
    fmindex_collection::occtable::naive::OccTable<256>
//    fmindex_collection::occtable::sdsl_wt_bldc::OccTable<256>,
//    fmindex_collection::occtable::sdsl_wt_epr::OccTable<256>
) {
    using OccTable = TestType;

   auto text = std::vector<uint8_t>{'H', 'a', 'l', 'l', 'o', ' ', 'W', 'e', 'l', 't',
                                    'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                    'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                    'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                    'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                    'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                    'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                    'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                    'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                    'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                    'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                    'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                    'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                    'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                    'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                    'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                    'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                    'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                    'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                    'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                    'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                    'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                    'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                    'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                    'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                    'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                    'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                    'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                    'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                    'x', 'y', 'z', 'x', 'y', 'z', 'x', 'y', 'z', 'x',
                                    };

    auto table = OccTable{text};

    REQUIRE(table.size() == text.size());

    SECTION("check that symbol() call works") {
        for (size_t i{0}; i < text.size(); ++i) {
            INFO(i);
            CHECK(table.symbol(i) == text.at(i));
        }
    }
    CHECK(!OccTable::name().empty());
    CHECK(!OccTable::extension().empty());
    CHECK(OccTable::Sigma == 256);

    SECTION("check all_ranks() is equal to prefix_rank() and rank()") {
        for (size_t idx{0}; idx < table.size(); ++idx) {
            auto [rank, prefix] = table.all_ranks(idx);
            for (size_t symb{1}; symb < 256; ++symb) {
                INFO(idx);
                INFO(symb);
                CHECK(rank[symb] == table.rank(idx, symb));
                CHECK(prefix[symb] == table.prefix_rank(idx, symb));
            }
        }
    }


}
