// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#pragma once

#include <fmindex-collection/bitvector/all.h>
#include <fmindex-collection/rankvector/rankvector.h>

#if SIZE_MAX == UINT64_MAX
#define AddIf64Bit(x) x,
#else
#define AddIf64Bit(x)
#endif
#define ALLRANKVECTORS_IMPL(Sigma) \
    fmindex_collection::rankvector::Naive<Sigma>, \
    (fmindex_collection::rankvector::MultiBitvector<Sigma, fmindex_collection::bitvector::Bitvector>), \
    (fmindex_collection::rankvector::MultiBitvector<Sigma, fmindex_collection::bitvector::CompactBitvector>), \
    (fmindex_collection::rankvector::MultiBitvector<Sigma, fmindex_collection::bitvector::CompactBitvector4Blocks>), \
    (fmindex_collection::rankvector::MultiBitvector<Sigma, fmindex_collection::bitvector::SparseBLEBitvector<>>), \
    (fmindex_collection::rankvector::MultiBitvector<Sigma, fmindex_collection::bitvector::SparseBLEBitvector<3>>), \
    (fmindex_collection::rankvector::MultiBitvector<Sigma, fmindex_collection::bitvector::SparseBLEBitvector<4>>), \
    (fmindex_collection::rankvector::MultiBitvector<Sigma, fmindex_collection::bitvector::SparseBLEBitvector<4, fmindex_collection::bitvector::Bitvector, fmindex_collection::bitvector::SparseBLEBitvector<>>>), \
    (fmindex_collection::rankvector::MultiBitvector<Sigma, fmindex_collection::bitvector::SparseBLEBitvector<5>>), \
    (fmindex_collection::rankvector::MultiBitvector<Sigma, fmindex_collection::bitvector::SparseBLEBitvector<6>>), \
    (fmindex_collection::rankvector::MultiBitvector<Sigma, fmindex_collection::bitvector::SparseBLEBitvector<6, fmindex_collection::bitvector::Bitvector, fmindex_collection::bitvector::SparseBLEBitvector<>>>), \
    (fmindex_collection::rankvector::MultiBitvector<Sigma, fmindex_collection::bitvector::SparseBLEBitvector<6, fmindex_collection::bitvector::Bitvector, fmindex_collection::bitvector::SparseBLEBitvector<3>>>), \
    (fmindex_collection::rankvector::MultiBitvector<Sigma, fmindex_collection::bitvector::SparseBLEBitvector<6, fmindex_collection::bitvector::Bitvector, fmindex_collection::bitvector::SparseBLEBitvector<4>>>), \
    (fmindex_collection::rankvector::MultiBitvector<Sigma, fmindex_collection::bitvector::SparseBLEBitvector<7>>), \
    (fmindex_collection::rankvector::MultiBitvector<Sigma, fmindex_collection::bitvector::SparseBLEBitvector<8>>), \
    fmindex_collection::rankvector::InterleavedBitvector8<Sigma>, \
    fmindex_collection::rankvector::InterleavedBitvector16<Sigma>, \
    fmindex_collection::rankvector::InterleavedBitvector32<Sigma>, \
    fmindex_collection::rankvector::InterleavedBitvector8Aligned<Sigma>, \
    fmindex_collection::rankvector::InterleavedBitvector16Aligned<Sigma>, \
    fmindex_collection::rankvector::InterleavedBitvector32Aligned<Sigma>, \
    fmindex_collection::rankvector::InterleavedEPR8<Sigma>, \
    fmindex_collection::rankvector::InterleavedEPR16<Sigma>, \
    fmindex_collection::rankvector::InterleavedEPR32<Sigma>, \
    fmindex_collection::rankvector::InterleavedEPR8Aligned<Sigma>, \
    fmindex_collection::rankvector::InterleavedEPR16Aligned<Sigma>, \
    fmindex_collection::rankvector::InterleavedEPR32Aligned<Sigma>, \
    fmindex_collection::rankvector::InterleavedEPRV2_8<Sigma>, \
    fmindex_collection::rankvector::InterleavedEPRV2_16<Sigma>, \
    AddIf64Bit(fmindex_collection::rankvector::InterleavedEPRV2_32<Sigma>) \
    fmindex_collection::rankvector::InterleavedEPRV2_8Aligned<Sigma>, \
    fmindex_collection::rankvector::InterleavedEPRV2_16Aligned<Sigma>, \
    AddIf64Bit(fmindex_collection::rankvector::InterleavedEPRV2_32Aligned<Sigma>) \
    fmindex_collection::rankvector::EPRV3_8<Sigma>, \
    fmindex_collection::rankvector::EPRV3_16<Sigma>, \
    AddIf64Bit(fmindex_collection::rankvector::EPRV3_32<Sigma>) \
    fmindex_collection::rankvector::EPRV4<Sigma>, \
    fmindex_collection::rankvector::EPRV5<Sigma>, \
    fmindex_collection::rankvector::DenseEPRV6<Sigma>, \
    fmindex_collection::rankvector::InterleavedEPRV7<Sigma>, \
    fmindex_collection::rankvector::InterleavedWavelet<Sigma>, \
    fmindex_collection::rankvector::Wavelet<Sigma>, \
    (fmindex_collection::rankvector::Wavelet<Sigma, fmindex_collection::bitvector::SparseBLEBitvector<>>), \
    fmindex_collection::rankvector::Double64ShortEPRV8<Sigma>, \
    fmindex_collection::rankvector::Double128ShortEPRV8<Sigma>, \
    fmindex_collection::rankvector::Double64EPRV8<Sigma>, \
    fmindex_collection::rankvector::Double128EPRV8<Sigma>, \
    fmindex_collection::rankvector::Double256EPRV8<Sigma>, \
    fmindex_collection::rankvector::Double512EPRV8<Sigma>
#if FMC_USE_SDSL
#define ALLRANKVECTORS(Sigma) \
    ALLRANKVECTORS_IMPL(Sigma), \
    fmindex_collection::rankvector::Sdsl_wt_bldc<Sigma>
#else
#define ALLRANKVECTORS(Sigma) ALLRANKVECTORS_IMPL(Sigma)
#endif

//!wt_epr is not working as expected
//    fmindex_collection::rankvector::Sdsl_wt_epr<255>,
