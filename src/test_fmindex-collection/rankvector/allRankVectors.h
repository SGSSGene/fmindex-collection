// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#pragma once

#include <fmindex-collection/bitvector/all.h>
#include <fmindex-collection/rankvector/rankvector.h>

#define ALLRANKVECTORS_IMPL \
    (fmindex_collection::rankvector::MultiBitvector<256, fmindex_collection::bitvector::SparseBLEBitvector<>>)

#if FMC_USE_SDSL
#define ALLRANKVECTORS \
    ALLRANKVECTORS_IMPL, \
    fmindex_collection::rankvector::Sdsl_wt_bldc<256>
#else
#define ALLRANKVECTORS ALLRANKVECTORS_IMPL
#endif

//!wt_epr is not working as expected
//    fmindex_collection::rankvector::Sdsl_wt_epr<255>,
