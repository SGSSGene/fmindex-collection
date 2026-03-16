// SPDX-FileCopyrightText: 2026 Simon Gene Gottlieb
// SPDX-License-Identifier: CC0-1.0
#include "../string/allStrings.h"
#include "../string/utils.h"

#include <fmindex-collection/fmindex/BiFMIndexKStep.h>
#include <fmindex-collection/search/SelectCursor.h>
#include <fmindex-collection/string/PairedFlattenedBitvectors2LPartialSymb.h>
#include <fmindex-collection/suffixarray/CSA.h>
#include <fstream>

/*template <size_t Sigma>
using String = fmc::string::PairedFlattenedBitvectors_512_64k<Sigma>;*/

template <size_t Sigma>
using String = fmc::string::PairedFlattenedBitvectorsPartialSymb_512_64k<Sigma>;

/*template <size_t Sigma>
using String = fmc::string::InterleavedBitvectorPrefix16<Sigma>;*/


TEST_CASE("benchmark fmindex extend left", "[fmindex][!benchmark][extend_left]") {

    static auto const& text = generateTexts<0, 4>();
    using Index      = fmc::BiFMIndex<4, String>::NoDelim;
    using IndexKStep = fmc::BiFMIndexKStep<4, String>::SetKStep<2>::NoDelim;
    static auto index       = Index{text, /*.samplingRate=*/1, /*.threadNbr=*/8};
    static auto index_kstep = IndexKStep{text, /*.samplingRate=*/1, /*.threadNbr=*/8};


    SECTION("benchmarking leftExtend") {
        auto rng = ankerl::nanobench::Rng{};

        auto bench = ankerl::nanobench::Bench{};
        bench.title("leftExtend")
             .relative(true)
             .batch(100);

        bench.run("BiFMIndex", [&]() {
            auto startI = rng.bounded(text[0].size()-101);
            auto cur = fmc::select_cursor_t<Index>{index};
            for (size_t i{startI+100}; i > startI; --i) {
                cur = cur.extendLeft(text[0][i]);
                assert(cur.count() > 0);
            }
            ankerl::nanobench::doNotOptimizeAway(cur);
        });
        bench.run("BiFMIndexKStep", [&]() {
            auto startI = rng.bounded(text[0].size()-101);
            auto cur = fmc::select_cursor_t<IndexKStep>{index_kstep};
            for (size_t i{startI+100}; i > startI; --i) {
                cur = cur.extendLeft(text[0][i]);
                assert(cur.count() > 0);
            }
            ankerl::nanobench::doNotOptimizeAway(cur);
        });

        bench.run("BiFMIndex-all", [&]() {
            auto startI = rng.bounded(text[0].size()-101);
            auto cur = fmc::select_cursor_t<Index>{index};
            for (size_t i{startI+100}; i > startI; --i) {
                auto cursors = cur.extendLeft();
                cur = cursors[text[0][i]];
                assert(cur.count() > 0);
            }
            ankerl::nanobench::doNotOptimizeAway(cur);
        });
        bench.run("BiFMIndex-all 2-step", [&]() {
            auto startI = rng.bounded(text[0].size()-101);
            auto cur = fmc::select_cursor_t<Index>{index};
            for (size_t i{startI+100}; i > startI; i -= 2) {
                auto cursors = cur.extendLeft();
                cur = cursors[text[0][i]];
                auto cursors2 = cur.extendLeft();
                cur = cursors2[text[0][i-1]];
                assert(cur.count() > 0);
            }
            ankerl::nanobench::doNotOptimizeAway(cur);
        });

        bench.run("BiFMIndexKStep-all", [&]() {
            auto startI = rng.bounded(text[0].size()-101);
            auto cur = fmc::select_cursor_t<IndexKStep>{index_kstep};
            for (size_t i{startI+100}; i > startI; --i) {
                auto cursors = cur.extendLeft();
                cur = cursors[text[0][i]];
                assert(cur.count() > 0);
            }
            ankerl::nanobench::doNotOptimizeAway(cur);
        });

        bench.run("BiFMIndexKStep-all 2-step", [&]() {
            auto startI = rng.bounded(text[0].size()-101);
            auto cur = fmc::select_cursor_t<IndexKStep>{index_kstep};
            for (size_t i{startI+100}; i > startI; i -= 2) {
                auto cursors = cur.extendLeftKStep();
                cur = cursors[text[0][i-1] + text[0][i]*4];
                assert(cur.count() > 0);
            }
            ankerl::nanobench::doNotOptimizeAway(cur);
        });


        bench.run("BiFMIndexKStep - 2step", [&]() {
            auto startI = rng.bounded(text[0].size()-101);
            auto cur = fmc::select_cursor_t<IndexKStep>{index_kstep};
            auto buf = std::array<size_t, 2>{};
            for (size_t i{startI+100}; i > startI; i -= 2) {
                buf = {text[0][i], text[0][i-1]};
                cur = cur.extendLeftKStep(buf);
            }
            ankerl::nanobench::doNotOptimizeAway(cur);
        });
        bench.run("BiFMIndex - 2step", [&]() {
            auto startI = rng.bounded(text[0].size()-101);
            auto cur = fmc::select_cursor_t<Index>{index};
            for (size_t i{startI+100}; i > startI; i -= 2) {
                cur = cur.extendLeft(text[0][i]);
                cur = cur.extendLeft(text[0][i-1]);
            }
            ankerl::nanobench::doNotOptimizeAway(cur);
        });

    }

    SECTION("benchmarking rightExtend") {
        auto rng = ankerl::nanobench::Rng{};

        auto bench = ankerl::nanobench::Bench{};
        bench.title("rightExtend")
             .relative(true)
             .batch(100);

        bench.run("BiFMIndex", [&]() {
            auto startI = rng.bounded(text[0].size()-101);
            auto cur = fmc::select_cursor_t<Index>{index};
            for (size_t i{startI}; i < startI+100; ++i) {
                cur = cur.extendRight(text[0][i]);
            }
            ankerl::nanobench::doNotOptimizeAway(cur);
        });
        bench.run("BiFMIndexKStep", [&]() {
            auto startI = rng.bounded(text[0].size()-101);
            auto cur = fmc::select_cursor_t<IndexKStep>{index_kstep};
            for (size_t i{startI}; i < startI+100; ++i) {
                cur = cur.extendRight(text[0][i]);
            }
            ankerl::nanobench::doNotOptimizeAway(cur);
        });

        bench.run("BiFMIndexKStep - 2step", [&]() {
            auto startI = rng.bounded(text[0].size()-101);
            auto cur = fmc::select_cursor_t<IndexKStep>{index_kstep};
            auto buf = std::array<size_t, 2>{};
            for (size_t i{startI}; i < startI+100; i += 2) {
                buf = {text[0][i], text[0][i+1]};
                cur = cur.extendRightKStep(buf);
            }
            ankerl::nanobench::doNotOptimizeAway(cur);
        });

    }

}
