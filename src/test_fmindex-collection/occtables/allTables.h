// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#pragma once

#include <fmindex-collection/occtable/all.h>

/*#define ALLTABLES \
%    fmindex_collection::occtable::bitvector::OccTable<256>, \
%    fmindex_collection::occtable::compactBitvector::OccTable<256>, \
    fmindex_collection::occtable::compactBitvectorPrefix::OccTable<256>, \
%    fmindex_collection::occtable::interleaved8::OccTable<256>, \
%    fmindex_collection::occtable::interleaved16::OccTable<256>, \
%    fmindex_collection::occtable::interleaved32::OccTable<256>, \
%    fmindex_collection::occtable::interleaved8Aligned::OccTable<256>, \
%    fmindex_collection::occtable::interleaved16Aligned::OccTable<256>, \
%    fmindex_collection::occtable::interleaved32Aligned::OccTable<256>, \
    fmindex_collection::occtable::interleavedPrefix::OccTable<256>, \
%    fmindex_collection::occtable::wavelet::OccTable<256>, \
%    fmindex_collection::occtable::interleavedWavelet::OccTable<256>, \
    fmindex_collection::occtable::interleavedWaveletAligned::OccTable<256>, \
    fmindex_collection::occtable::interleavedWavelet32::OccTable<256>, \
    fmindex_collection::occtable::interleavedWavelet32Aligned::OccTable<256>, \
%    fmindex_collection::occtable::interleavedEPR8::OccTable<256>, \
%    fmindex_collection::occtable::interleavedEPR16::OccTable<256>, \
%    fmindex_collection::occtable::interleavedEPR32::OccTable<256>, \
%    fmindex_collection::occtable::interleavedEPR8Aligned::OccTable<256>, \
%    fmindex_collection::occtable::interleavedEPR16Aligned::OccTable<256>, \
%    fmindex_collection::occtable::interleavedEPR32Aligned::OccTable<256>, \
%    fmindex_collection::occtable::interleavedEPR8V2::OccTable<256>, \
%    fmindex_collection::occtable::interleavedEPR16V2::OccTable<256>, \
%    fmindex_collection::occtable::interleavedEPR32V2::OccTable<256>, \
%    fmindex_collection::occtable::epr8V3::OccTable<256>, \
%    fmindex_collection::occtable::epr16V3::OccTable<256>, \
%    fmindex_collection::occtable::epr32V3::OccTable<256>, \
%    fmindex_collection::occtable::eprV4::OccTable<256>, \
%    fmindex_collection::occtable::eprV5::OccTable<256>, \
%    fmindex_collection::occtable::eprV6::OccTable<256>, \
%    fmindex_collection::occtable::rlebwt::OccTable<256>, \
%    fmindex_collection::occtable::naive::OccTable<256>
%    fmindex_collection::occtable::sdsl_wt_bldc::OccTable<256>
    fmindex_collection::occtable::sdsl_wt_epr::OccTable<256>
*/

#if UINT64_MAX == SIZE_MAX
    #define setif64bit(V) V,
#else
    #define setif64bit(V)
#endif

#define ALLTABLES_IMPL \
    fmindex_collection::occtable::Naive<256>, \
    fmindex_collection::occtable::Bitvector<256>, \
    fmindex_collection::occtable::L1Bitvector<256>, \
    fmindex_collection::occtable::CompactBitvector<256>, \
    fmindex_collection::occtable::CompactBitvector4Blocks<256>, \
    fmindex_collection::occtable::Interleaved_8<256>, \
    fmindex_collection::occtable::Interleaved_16<256>, \
    setif64bit(fmindex_collection::occtable::Interleaved_32<256>) \
    fmindex_collection::occtable::Interleaved_8Aligned<256>, \
    fmindex_collection::occtable::Interleaved_16Aligned<256>, \
    setif64bit(fmindex_collection::occtable::Interleaved_32Aligned<256>) \
    fmindex_collection::occtable::Epr_8<256>, \
    fmindex_collection::occtable::Epr_16<256>, \
    setif64bit(fmindex_collection::occtable::Epr_32<256>) \
    fmindex_collection::occtable::Epr_8Aligned<256>, \
    fmindex_collection::occtable::Epr_16Aligned<256>, \
    setif64bit(fmindex_collection::occtable::Epr_32Aligned<256>) \
    fmindex_collection::occtable::EprV2_8<256>, \
    fmindex_collection::occtable::EprV2_16<256>, \
    setif64bit(fmindex_collection::occtable::EprV2_32<256>) \
    fmindex_collection::occtable::EprV2_8Aligned<256>, \
    fmindex_collection::occtable::EprV2_16Aligned<256>, \
    setif64bit(fmindex_collection::occtable::EprV2_32Aligned<256>) \
    fmindex_collection::occtable::EprV3_8<256>, \
    fmindex_collection::occtable::EprV3_16<256>, \
    setif64bit(fmindex_collection::occtable::EprV3_32<256>) \
    fmindex_collection::occtable::EprV4<256>, \
    fmindex_collection::occtable::EprV5<256>, \
    fmindex_collection::occtable::EprV6<256>, \
    fmindex_collection::occtable::EprV7<256>, \
    fmindex_collection::occtable::InterleavedWavelet<256>, \
    fmindex_collection::occtable::Wavelet<256>, \
    fmindex_collection::occtable::RuntimeLengthEncoded2<256>, \
    fmindex_collection::occtable::RuntimeLengthEncoded3<256>, \
    fmindex_collection::occtable::RuntimeLengthEncoded4<256>, \
    fmindex_collection::occtable::RecursiveRuntimeLengthEncodedD2<256>

#if FMC_USE_SDSL
#define ALLTABLES \
    ALLTABLES_IMPL, \
    fmindex_collection::occtable::Sdsl_wt_bldc<256>
#else
#define ALLTABLES ALLTABLES_IMPL
#endif
