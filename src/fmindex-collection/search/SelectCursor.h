// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../fmindex/BiFMIndexCursor.h"
#include "../fmindex/FMIndexCursor.h"
#include "../fmindex/RBiFMIndexCursor.h"
#include "../fmindex/ReverseFMIndexCursor.h"


namespace fmindex_collection {

template <typename Index>
struct SelectIndexCursor;

template <typename OccTable, typename TCSA>
struct SelectIndexCursor<BiFMIndex<OccTable, TCSA>> {
    using cursor_t = BiFMIndexCursor<BiFMIndex<OccTable, TCSA>>;
};

template <typename OccTable, typename TCSA>
struct SelectIndexCursor<FMIndex<OccTable, TCSA>> {
    using cursor_t = FMIndexCursor<FMIndex<OccTable, TCSA>>;
};

template <typename OccTable, typename TCSA>
struct SelectIndexCursor<RBiFMIndex<OccTable, TCSA>> {
    using cursor_t = RBiFMIndexCursor<RBiFMIndex<OccTable, TCSA>>;
};

template <typename OccTable, typename TCSA>
struct SelectIndexCursor<ReverseFMIndex<OccTable, TCSA>> {
    using cursor_t = ReverseFMIndexCursor<ReverseFMIndex<OccTable, TCSA>>;
};


template <typename Index>
struct SelectLeftIndexCursor;

template <typename OccTable, typename TCSA>
struct SelectLeftIndexCursor<BiFMIndex<OccTable, TCSA>> {
    using cursor_t = LeftBiFMIndexCursor<BiFMIndex<OccTable, TCSA>>;
};

template <typename OccTable, typename TCSA>
struct SelectLeftIndexCursor<FMIndex<OccTable, TCSA>> {
    using cursor_t = FMIndexCursor<FMIndex<OccTable, TCSA>>;
};

template <typename OccTable, typename TCSA>
struct SelectLeftIndexCursor<RBiFMIndex<OccTable, TCSA>> {
    using cursor_t = LeftRBiFMIndexCursor<RBiFMIndex<OccTable, TCSA>>;
};


template <typename Index>
using select_cursor_t      = typename SelectIndexCursor<Index>::cursor_t;

template <typename Index>
using select_left_cursor_t = typename SelectLeftIndexCursor<Index>::cursor_t;

}
