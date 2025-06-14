// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../fmindex/BiFMIndexCursor.h"
#include "../fmindex/FMIndexCursor.h"
#include "../fmindex/KMerFMIndexCursor.h"
#include "../fmindex/MirroredBiFMIndexCursor.h"
#include "../fmindex/ReverseFMIndexCursor.h"


namespace fmindex_collection {

template <typename Index>
struct SelectIndexCursor;

template <typename String, typename TCSA>
struct SelectIndexCursor<BiFMIndex<String, TCSA>> {
    using cursor_t = BiFMIndexCursor<BiFMIndex<String, TCSA>>;
};

template <typename String, typename TCSA>
struct SelectIndexCursor<FMIndex<String, TCSA>> {
    using cursor_t = FMIndexCursor<FMIndex<String, TCSA>>;
};

template <typename String, size_t KMer, typename TCSA>
struct SelectIndexCursor<KMerFMIndex<String, KMer, TCSA>> {
    using cursor_t = KMerFMIndexCursor<KMerFMIndex<String, KMer, TCSA>>;
};

template <typename String, typename TCSA>
struct SelectIndexCursor<MirroredBiFMIndex<String, TCSA>> {
    using cursor_t = MirroredBiFMIndexCursor<MirroredBiFMIndex<String, TCSA>>;
};

template <typename String, typename TCSA>
struct SelectIndexCursor<ReverseFMIndex<String, TCSA>> {
    using cursor_t = ReverseFMIndexCursor<ReverseFMIndex<String, TCSA>>;
};


template <typename Index>
struct SelectLeftIndexCursor;

template <typename String, typename TCSA>
struct SelectLeftIndexCursor<BiFMIndex<String, TCSA>> {
    using cursor_t = LeftBiFMIndexCursor<BiFMIndex<String, TCSA>>;
};

template <typename String, typename TCSA>
struct SelectLeftIndexCursor<FMIndex<String, TCSA>> {
    using cursor_t = FMIndexCursor<FMIndex<String, TCSA>>;
};

template <typename String, size_t KMer, typename TCSA>
struct SelectLeftIndexCursor<KMerFMIndex<String, KMer, TCSA>> {
    using cursor_t = KMerFMIndexCursor<KMerFMIndex<String, KMer, TCSA>>;
};

template <typename String, typename TCSA>
struct SelectLeftIndexCursor<MirroredBiFMIndex<String, TCSA>> {
    using cursor_t = LeftMirroredBiFMIndexCursor<MirroredBiFMIndex<String, TCSA>>;
};


template <typename Index>
using select_cursor_t      = typename SelectIndexCursor<Index>::cursor_t;

template <typename Index>
using select_left_cursor_t = typename SelectLeftIndexCursor<Index>::cursor_t;

}
