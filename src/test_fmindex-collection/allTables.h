// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include <fmindex-collection/occtable/all.h>

#define ALLTABLES \
    fmindex_collection::occtable::bitvector::OccTable<256>, \
    fmindex_collection::occtable::compactBitvector::OccTable<256>, \
    fmindex_collection::occtable::compactBitvectorPrefix::OccTable<256>, \
    fmindex_collection::occtable::interleaved8::OccTable<256>, \
    fmindex_collection::occtable::interleaved16::OccTable<256>, \
    fmindex_collection::occtable::interleaved32::OccTable<256>, \
    fmindex_collection::occtable::interleaved8Aligned::OccTable<256>, \
    fmindex_collection::occtable::interleaved16Aligned::OccTable<256>, \
    fmindex_collection::occtable::interleaved32Aligned::OccTable<256>, \
    fmindex_collection::occtable::interleavedPrefix::OccTable<256>, \
    fmindex_collection::occtable::wavelet::OccTable<256>, \
    fmindex_collection::occtable::interleavedWavelet::OccTable<256>, \
    fmindex_collection::occtable::interleavedWaveletAligned::OccTable<256>, \
    fmindex_collection::occtable::interleavedWavelet32::OccTable<256>, \
    fmindex_collection::occtable::interleavedWavelet32Aligned::OccTable<256>, \
    fmindex_collection::occtable::interleavedEPR8::OccTable<256>, \
    fmindex_collection::occtable::interleavedEPR16::OccTable<256>, \
    fmindex_collection::occtable::interleavedEPR32::OccTable<256>, \
    fmindex_collection::occtable::interleavedEPR8Aligned::OccTable<256>, \
    fmindex_collection::occtable::interleavedEPR16Aligned::OccTable<256>, \
    fmindex_collection::occtable::interleavedEPR32Aligned::OccTable<256>, \
    fmindex_collection::occtable::interleavedEPR8V2::OccTable<256>, \
    fmindex_collection::occtable::interleavedEPR16V2::OccTable<256>, \
    fmindex_collection::occtable::interleavedEPR32V2::OccTable<256>, \
    fmindex_collection::occtable::epr8V3::OccTable<256>, \
    fmindex_collection::occtable::epr16V3::OccTable<256>, \
    fmindex_collection::occtable::epr32V3::OccTable<256>, \
    fmindex_collection::occtable::eprV4::OccTable<256>, \
    fmindex_collection::occtable::eprV5::OccTable<256>, \
    fmindex_collection::occtable::eprV6::OccTable<256>, \
    fmindex_collection::occtable::naive::OccTable<256>
//    fmindex_collection::occtable::sdsl_wt_bldc::OccTable<256>
//    fmindex_collection::occtable::sdsl_wt_epr::OccTable<256>

