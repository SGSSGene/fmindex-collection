// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#pragma once

#include <fmindex-collection/bitvector/all.h>
#include <fmindex-collection/rankvector/rankvector.h>

#define ALLSYMBOLVECTORS_IMPL \
    fmindex_collection::rankvector::Naive<256>, \
    (fmindex_collection::rankvector::MultiBitvector<256, fmindex_collection::bitvector::Bitvector>), \
    (fmindex_collection::rankvector::MultiBitvector<256, fmindex_collection::bitvector::CompactBitvector>), \
    (fmindex_collection::rankvector::MultiBitvector<256, fmindex_collection::bitvector::CompactBitvector4Blocks>), \
    fmindex_collection::rankvector::InterleavedBitvector8<256>, \
    fmindex_collection::rankvector::InterleavedBitvector16<256>, \
    fmindex_collection::rankvector::InterleavedBitvector32<256>, \
    fmindex_collection::rankvector::InterleavedBitvector8Aligned<256>, \
    fmindex_collection::rankvector::InterleavedBitvector16Aligned<256>, \
    fmindex_collection::rankvector::InterleavedBitvector32Aligned<256>, \
    fmindex_collection::rankvector::InterleavedEPR8<256>, \
    fmindex_collection::rankvector::InterleavedEPR16<256>, \
    fmindex_collection::rankvector::InterleavedEPR32<256>, \
    fmindex_collection::rankvector::InterleavedEPR8Aligned<256>, \
    fmindex_collection::rankvector::InterleavedEPR16Aligned<256>, \
    fmindex_collection::rankvector::InterleavedEPR32Aligned<256>, \
    fmindex_collection::rankvector::InterleavedEPRV2_8<256>, \
    fmindex_collection::rankvector::InterleavedEPRV2_16<256>, \
    fmindex_collection::rankvector::InterleavedEPRV2_32<256>, \
    fmindex_collection::rankvector::InterleavedEPRV2_8Aligned<256>, \
    fmindex_collection::rankvector::InterleavedEPRV2_16Aligned<256>, \
    fmindex_collection::rankvector::InterleavedEPRV2_32Aligned<256>, \
    fmindex_collection::rankvector::EPRV3_8<256>, \
    fmindex_collection::rankvector::EPRV3_16<256>, \
    fmindex_collection::rankvector::EPRV3_32<256>, \
    fmindex_collection::rankvector::EPRV4<256>, \
    fmindex_collection::rankvector::EPRV5<256>, \
    fmindex_collection::rankvector::DenseEPRV6<256>, \
    fmindex_collection::rankvector::InterleavedEPRV7<256>, \
    fmindex_collection::rankvector::InterleavedWavelet<256>, \
    fmindex_collection::rankvector::Wavelet<256>, \
    fmindex_collection::rankvector::RLEInstance<256>

#if FMC_USE_SDSL
#define ALLSYMBOLVECTORS \
    ALLSYMBOLVECTORS_IMPL, \
    fmindex_collection::rankvector::Sdsl_wt_bldc<256>
#else
#define ALLSYMBOLVECTORS ALLSYMBOLVECTORS_IMPL

#endif


//!wt_epr is not working as expected
//    fmindex_collection::rankvector::Sdsl_wt_epr<255>,


