#include <fmindex-collection/occtable/all.h>
#include <fmindex-collection/BiFMIndex.h>
#include <fmindex-collection/search/SearchPseudo.h>
#include <search_schemes/generator/all.h>

#include <catch2/catch.hpp>

TEMPLATE_TEST_CASE("searching with PseudoSearch", "[search]",
    fmindex_collection::occtable::bitvector::OccTable<256>,
    fmindex_collection::occtable::compactBitvector::OccTable<256>,
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

    auto input  = std::vector<uint8_t>{'A', 'A', 'A', 'C', 'A', 'A', 'A', 'C', 'A', 'A', 'A'};

    auto index = fmindex_collection::BiFMIndex<OccTable>{std::vector<std::vector<uint8_t>>{input}, 1};

    SECTION("check symbol call to occurrence table") {
        REQUIRE(input.size()+1 == index.size());
        CHECK(index.occ.symbol( 0) == 'A');
        CHECK(index.occ.symbol( 1) == 'A');
        CHECK(index.occ.symbol( 2) == 'A');
        CHECK(index.occ.symbol( 3) == 'C');
        CHECK(index.occ.symbol( 4) == 'C');
        CHECK(index.occ.symbol( 5) == '\0');
        CHECK(index.occ.symbol( 6) == 'A');
        CHECK(index.occ.symbol( 7) == 'A');
        CHECK(index.occ.symbol( 8) == 'A');
        CHECK(index.occ.symbol( 9) == 'A');
        CHECK(index.occ.symbol(10) == 'A');
        CHECK(index.occ.symbol(11) == 'A');
    }

    auto query = std::vector<std::vector<uint8_t>> {std::vector<uint8_t>{'A'}};
    auto search_scheme = search_schemes::generator::backtracking(1, 0, 0);
    fmindex_collection::search_pseudo::search<true>(index, query, search_scheme, [](auto qidx, auto result, auto errors) {
        CHECK(qidx == 0);
        CHECK(errors == 0);
        CHECK(result.lb == 1);
        CHECK(result.count() == 9);
    });

}
TEMPLATE_TEST_CASE("searching with collection and PseudoSearch", "[collection]",
    fmindex_collection::occtable::bitvector::OccTable<256>,
    fmindex_collection::occtable::compactBitvector::OccTable<256>,
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

    auto input  = std::vector<std::vector<uint8_t>>{{'A', 'A', 'A', 'C', 'A', 'A', 'A', 'C', 'A', 'A', 'A'},
                                                    {'A', 'A', 'A', 'B', 'A', 'A', 'A', 'B', 'A', 'A', 'A'}};

    auto index = fmindex_collection::BiFMIndex<OccTable>{input, 1};

    SECTION("check symbol call to occurrence table") {
        auto expected = std::vector<uint8_t>{'A', 'A', 'A', 'A', 'A', 'A', 'B', 'C', 'B', '\0', 'C', '\0',
                                             'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A', 'A'};
        REQUIRE(index.size() == expected.size());
        for (size_t i{0}; i < expected.size(); ++i) {
            INFO(i);
            CHECK(index.occ.symbol(i) == expected[i]);
        }
    }

    auto query = std::vector<std::vector<uint8_t>> {std::vector<uint8_t>{'A'}};
    auto search_scheme = search_schemes::generator::backtracking(1, 0, 0);
    fmindex_collection::search_pseudo::search<true>(index, query, search_scheme, [](auto qidx, auto result, auto errors) {
        CHECK(qidx == 0);
        CHECK(errors == 0);
        CHECK(result.lb == 2);
        CHECK(result.count() == 18);
    });

    auto expected = std::vector<std::tuple<size_t, size_t>> {
        std::make_tuple(1ul, 11ul),
        std::make_tuple(0ul, 11ul),
        std::make_tuple(1ul, 10ul),
        std::make_tuple(0ul, 10ul),
        std::make_tuple(1ul, 9ul),
        std::make_tuple(0ul, 9ul),
        std::make_tuple(1ul, 8ul),
        std::make_tuple(0ul, 8ul),
        std::make_tuple(1ul, 4ul),
        std::make_tuple(1ul, 0ul),
        std::make_tuple(0ul, 4ul),
        std::make_tuple(0ul, 0ul),
        std::make_tuple(1ul, 5ul),
        std::make_tuple(1ul, 1ul),
        std::make_tuple(0ul, 5ul),
        std::make_tuple(0ul, 1ul),
        std::make_tuple(1ul, 6ul),
        std::make_tuple(1ul, 2ul),
        std::make_tuple(0ul, 6ul),
        std::make_tuple(0ul, 2ul),
        std::make_tuple(1ul, 7ul),
        std::make_tuple(1ul, 3ul),
        std::make_tuple(0ul, 7ul),
        std::make_tuple(0ul, 3ul)
    };

    for (size_t i{0}; i < expected.size(); ++i) {
        INFO(i);
        auto [il, pl] = index.locate(i);
        auto [ir, pr] = expected[i];
        CHECK(il == ir);
        CHECK(pl == pr);
    }
}
