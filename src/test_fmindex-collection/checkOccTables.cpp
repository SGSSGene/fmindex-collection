
#include <fmindex-collection/occtable/all.h>
#include <catch2/catch.hpp>


TEMPLATE_TEST_CASE("check if occ table is working", "[OccTable]",
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

    auto text = std::vector<uint8_t>{'H', 'a', 'l', 'l', 'o', ' ', 'W', 'e', 'l', 't'};
    auto table = OccTable{text};

    SECTION("check that symbol() call works") {
        REQUIRE(table.size() == text.size());
        for (size_t i{0}; i < text.size(); ++i) {
            CHECK(table.symbol(i) == text.at(i));
        }
    }
    SECTION("test complete table 'H'") {
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
}
