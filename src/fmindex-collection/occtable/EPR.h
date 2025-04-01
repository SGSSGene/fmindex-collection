// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../string/all.h"
#include "../bitvector/all.h"
#include "GenericOccTable.h"

namespace fmindex_collection::occtable {
/*!TODO Still needs a few more occtables translated:
 *
 *  Nice to have, but removed at some point
 *  - compactBitvectorPrefix
 *  - InterleavedPrefix
 *  - sdsl_wt_epr
 *  - InterleavedWaveletAligned
 *  - InterleavedWavelet32 (-/Aligned)
 */

template <size_t TSigma>
using Naive = GenericOccTable<string::Naive<TSigma>, "Naive", "n">;
static_assert(checkOccTable<Naive>);

template <size_t TSigma>
using Bitvector = GenericOccTable<string::MultiBitvector<TSigma, fmindex_collection::bitvector::Bitvector>, "bitvector", "bv">;
static_assert(checkOccTable<Bitvector>);

template <size_t TSigma>
using L1Bitvector = GenericOccTable<string::MultiBitvector<TSigma, fmindex_collection::bitvector::L1Bitvector>, "l1bitvector", "l1bv">;
static_assert(checkOccTable<L1Bitvector>);

template <size_t TSigma>
using CompactBitvector = GenericOccTable<string::MultiBitvector<TSigma, fmindex_collection::bitvector::CompactBitvector>, "compact bitvector", "cbv">;
static_assert(checkOccTable<CompactBitvector>);

template <size_t TSigma>
using CompactBitvector4Blocks = GenericOccTable<string::MultiBitvector<TSigma, fmindex_collection::bitvector::CompactBitvector4Blocks>, "compact bitvector 4 blocks", "cbv4">;
static_assert(checkOccTable<CompactBitvector4Blocks>);

template <size_t TSigma>
using Interleaved_8 = GenericOccTable<string::InterleavedBitvector8<TSigma>, "interleaved 8", "i8">;
static_assert(checkOccTable<Interleaved_8>);

template <size_t TSigma>
using Interleaved_16 = GenericOccTable<string::InterleavedBitvector16<TSigma>, "interleaved 16", "i16">;
static_assert(checkOccTable<Interleaved_16>);

#if SIZE_MAX == UINT64_MAX
template <size_t TSigma>
using Interleaved_32 = GenericOccTable<string::InterleavedBitvector32<TSigma>, "interleaved 32", "i32">;
static_assert(checkOccTable<Interleaved_32>);
#endif

template <size_t TSigma>
using Interleaved_8Aligned = GenericOccTable<string::InterleavedBitvector8Aligned<TSigma>, "interleaved 8", "i8">;
static_assert(checkOccTable<Interleaved_8Aligned>);

template <size_t TSigma>
using Interleaved_16Aligned = GenericOccTable<string::InterleavedBitvector16Aligned<TSigma>, "interleaved 16", "i16">;
static_assert(checkOccTable<Interleaved_16Aligned>);

#if SIZE_MAX == UINT64_MAX
template <size_t TSigma>
using Interleaved_32Aligned = GenericOccTable<string::InterleavedBitvector32Aligned<TSigma>, "interleaved 32", "i32">;
static_assert(checkOccTable<Interleaved_32Aligned>);
#endif

template <size_t TSigma>
using Epr_8 = GenericOccTable<string::InterleavedEPR8<TSigma>, "Epr (8bit)", "epr_8">;
static_assert(checkOccTable<Epr_8>);

template <size_t TSigma>
using Epr_16 = GenericOccTable<string::InterleavedEPR16<TSigma>, "Epr (16bit)", "epr_16">;
static_assert(checkOccTable<Epr_16>);

#if SIZE_MAX == UINT64_MAX
template <size_t TSigma>
using Epr_32 = GenericOccTable<string::InterleavedEPR32<TSigma>, "Epr (32bit)", "epr_32">;
static_assert(checkOccTable<Epr_32>);
#endif

template <size_t TSigma>
using Epr_8Aligned = GenericOccTable<string::InterleavedEPR8Aligned<TSigma>, "Epr (8bit, aligned)", "epr_8a">;
static_assert(checkOccTable<Epr_8Aligned>);

template <size_t TSigma>
using Epr_16Aligned = GenericOccTable<string::InterleavedEPR16Aligned<TSigma>, "Epr (16bit, aligned)", "epr_16a">;
static_assert(checkOccTable<Epr_16Aligned>);

#if SIZE_MAX == UINT64_MAX
template <size_t TSigma>
using Epr_32Aligned = GenericOccTable<string::InterleavedEPR32Aligned<TSigma>, "Epr (32bit, aligned)", "epr_32a">;
static_assert(checkOccTable<Epr_32Aligned>);
#endif

template <size_t TSigma>
using EprV2_8 = GenericOccTable<string::InterleavedEPRV2_8<TSigma>, "EprV2 (8bit)", "eprv2_8">;
static_assert(checkOccTable<EprV2_8>);

template <size_t TSigma>
using EprV2_16 = GenericOccTable<string::InterleavedEPRV2_16<TSigma>, "EprV2 (16bit)", "eprv2_16">;
static_assert(checkOccTable<EprV2_16>);

#if SIZE_MAX == UINT64_MAX
template <size_t TSigma>
using EprV2_32 = GenericOccTable<string::InterleavedEPRV2_32<TSigma>, "EprV2 (32bit)", "eprv2_32">;
static_assert(checkOccTable<EprV2_32>);
#endif

template <size_t TSigma>
using EprV2_8Aligned = GenericOccTable<string::InterleavedEPRV2_8Aligned<TSigma>, "EprV2 (8bit, aligned)", "eprv2_8a">;
static_assert(checkOccTable<EprV2_8Aligned>);

template <size_t TSigma>
using EprV2_16Aligned = GenericOccTable<string::InterleavedEPRV2_16Aligned<TSigma>, "EprV2 (16bit, aligned)", "eprv2_16a">;
static_assert(checkOccTable<EprV2_16Aligned>);

#if SIZE_MAX == UINT64_MAX
template <size_t TSigma>
using EprV2_32Aligned = GenericOccTable<string::InterleavedEPRV2_32Aligned<TSigma>, "EprV2 (32bit, aligned)", "eprv2_32a">;
static_assert(checkOccTable<EprV2_32Aligned>);
#endif

template <size_t TSigma>
using EprV3_8 = GenericOccTable<string::EPRV3_8<TSigma>, "EPRV3 (8bit)", "eprv3_8">;
static_assert(checkOccTable<EprV3_8>);

template <size_t TSigma>
using EprV3_16 = GenericOccTable<string::EPRV3_16<TSigma>, "EPRV3 (16bit)", "eprv3_16">;
static_assert(checkOccTable<EprV3_16>);

#if SIZE_MAX == UINT64_MAX
template <size_t TSigma>
using EprV3_32 = GenericOccTable<string::EPRV3_32<TSigma>, "EPRV3 (32bit)", "eprv3_32">;
static_assert(checkOccTable<EprV3_32>);
#endif

template <size_t TSigma>
using EprV4 = GenericOccTable<string::EPRV4<TSigma>, "EPRV4", "eprv4">;
static_assert(checkOccTable<EprV4>);

template <size_t TSigma>
using EprV5 = GenericOccTable<string::EPRV5<TSigma>, "EPRV5", "eprv5">;
static_assert(checkOccTable<EprV5>);

template <size_t TSigma>
using EprV6 = GenericOccTable<string::DenseEPRV6<TSigma>, "EPRV6", "eprv6">;
static_assert(checkOccTable<EprV6>);

template <size_t TSigma>
using EprV7 = GenericOccTable<string::InterleavedEPRV7<TSigma>, "EPRV7", "eprv7">;
static_assert(checkOccTable<EprV7>);

template <size_t TSigma>
using Wavelet = GenericOccTable<string::Wavelet<TSigma>, "Wavelet", "w">;
static_assert(checkOccTable<Wavelet>);

template <size_t TSigma>
using InterleavedWavelet = GenericOccTable<string::InterleavedWavelet<TSigma>, "Interleaved Wavelet", "iw">;
static_assert(checkOccTable<InterleavedWavelet>);

/**
 * A runlength encoded bwt, see unpublished paper
 *
 */
template <size_t TSigma>
using RunBlockEncoded2 = GenericOccTable<string::RBBwt<TSigma, 2>, "run block encoded 2", "rbbwt">;
static_assert(checkOccTable<RunBlockEncoded2>);

template <size_t TSigma>
using RunBlockEncoded3 = GenericOccTable<string::RBBwt<TSigma, 3>, "run block encoded 3", "rbbwt3">;
static_assert(checkOccTable<RunBlockEncoded3>);

template <size_t TSigma>
using RunBlockEncoded4 = GenericOccTable<string::RBBwt<TSigma, 4>, "run block encoded 4", "rbbwt4">;
static_assert(checkOccTable<RunBlockEncoded4>);

/**
 * A run length encoded bwt, see unpublished paper, but with recursive bitvectors
 *
 */
template <size_t TSigma>
using RecursiveRunBlockEncodedD2 = GenericOccTable<string::rRBBwt<TSigma, 2>, "recursive run block encoded, depth 2", "rec_rbbwt2">;
static_assert(checkOccTable<RecursiveRunBlockEncodedD2>);

#ifdef FMC_USE_SDSL

template <size_t TSigma>
using Sdsl_wt_bldc = GenericOccTable<string::Sdsl_wt_bldc<TSigma>, "SDSL wt_bldc", "sdsl_wt_bldc">;
static_assert(checkOccTable<Sdsl_wt_bldc>);

#endif
}
