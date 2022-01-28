#include <fmindex-collection/occtable/all.h>
#include <fmindex-collection/BiFMIndex.h>
#include <fmindex-collection/search/SearchPseudo.h>
#include <search_schemes/generator/all.h>

#include <catch2/catch.hpp>
#include <divsufsort64.h>

TEMPLATE_TEST_CASE("searching with PseudoSearch", "[search]",
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

    auto input  = std::vector<uint8_t>{'A', 'A', 'A', 'C', 'A', 'A', 'A', 'C', 'A', 'A', 'A'};

    auto index = fmindex_collection::BiFMIndex<OccTable>{input, 1};

    SECTION("check symbol call to occurence table") {
        REQUIRE(input.size() == index.size());
        CHECK(index.occ.symbol( 0) == 'A');
        CHECK(index.occ.symbol( 1) == 'A');
        CHECK(index.occ.symbol( 2) == 'C');
        CHECK(index.occ.symbol( 3) == 'C');
        CHECK(index.occ.symbol( 4) == 'A');
        CHECK(index.occ.symbol( 5) == 'A');
        CHECK(index.occ.symbol( 6) == 'A');
        CHECK(index.occ.symbol( 7) == 'A');
        CHECK(index.occ.symbol( 8) == 'A');
        CHECK(index.occ.symbol( 9) == 'A');
        CHECK(index.occ.symbol(10) == 'A');
    }

    auto query = std::vector<std::vector<uint8_t>> {std::vector<uint8_t>{'A'}};
    auto search_scheme = search_schemes::generator::backtracking(1, 0, 0);
    fmindex_collection::search_pseudo::search<true>(index, query, search_scheme, [](auto qidx, auto result, auto errors) {
        CHECK(qidx == 0);
        CHECK(errors == 0);
        CHECK(result.lb == 0);
        CHECK(result.count() == 9);
    });


/*    auto sa = std::vector<int64_t>{};
    sa.resize(input.size());
    auto error = divsufsort64(static_cast<uint8_t const*>(input.data()), sa.data(), input.size());
    if (error != 0) {
        throw std::runtime_error("some error while creating the suffix array");
    }

    CHECK(sa[ 0] == 10);
    CHECK(sa[ 1] ==  9);
    CHECK(sa[ 2] ==  8);
    CHECK(sa[ 3] ==  4);
    CHECK(sa[ 4] ==  0);
    CHECK(sa[ 5] ==  5);
    CHECK(sa[ 6] ==  1);
    CHECK(sa[ 7] ==  6);
    CHECK(sa[ 8] ==  2);
    CHECK(sa[ 9] ==  7);
    CHECK(sa[10] ==  3);



    auto bwt    = std::vector<uint8_t>{'t', 'o', '\0', ' ', 'H', 'W', 'a', 'l', 'e', 'l', 'l'};
    auto bwtRev = std::vector<uint8_t>{'H', 'W', 'a', 'e', 'l', 'l', 'l', 't', 'o', ' ', '\0'};
    auto sa     = std::vector<size_t>{ 10, 5, 0,  6,  1,  7,  2,  3,  8,  4,  9 };*/

}
