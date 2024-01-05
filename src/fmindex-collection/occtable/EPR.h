// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "../rankvector/rankvector.h"
#include "../bitvector/bitvector.h"
#include "GenericOccTable.h"

namespace fmindex_collection::occtable {
/*!TODO Still needs a few more occtables translated:
 *  - rlebwt
 *
 *  - compactBitvectorPrefix
 *  - InterleavedPrefix
 *  - sdsl_wt_epr
 *  - InterleavedWaveletAligned
 *  - InterleavedWavelet32 (-/Aligned)
 */

template <uint64_t TSigma>
using Naive = GenericOccTable<rankvector::Naive<TSigma>, "Naive", "n">;
static_assert(checkOccTable<Naive>);

template <uint64_t TSigma>
using Bitvector = GenericOccTable<rankvector::MultiBitvector<TSigma, fmindex_collection::bitvector::Bitvector>, "bitvector", "bv">;
static_assert(checkOccTable<Bitvector>);

template <uint64_t TSigma>
using CompactBitvector = GenericOccTable<rankvector::MultiBitvector<TSigma, fmindex_collection::bitvector::CompactBitvector>, "compact bitvector", "cbv">;
static_assert(checkOccTable<CompactBitvector>);

template <uint64_t TSigma>
using CompactBitvector4Blocks = GenericOccTable<rankvector::MultiBitvector<TSigma, fmindex_collection::bitvector::CompactBitvector4Blocks>, "compact bitvector 4 blocks", "cbv4">;
static_assert(checkOccTable<CompactBitvector4Blocks>);

template <uint64_t TSigma>
using Interleaved_8 = GenericOccTable<rankvector::InterleavedBitvector8<TSigma>, "interleaved 8", "i8">;
static_assert(checkOccTable<Interleaved_8>);

template <uint64_t TSigma>
using Interleaved_16 = GenericOccTable<rankvector::InterleavedBitvector16<TSigma>, "interleaved 16", "i16">;
static_assert(checkOccTable<Interleaved_16>);

template <uint64_t TSigma>
using Interleaved_32 = GenericOccTable<rankvector::InterleavedBitvector32<TSigma>, "interleaved 32", "i32">;
static_assert(checkOccTable<Interleaved_32>);

template <uint64_t TSigma>
using Interleaved_8Aligned = GenericOccTable<rankvector::InterleavedBitvector8Aligned<TSigma>, "interleaved 8", "i8">;
static_assert(checkOccTable<Interleaved_8Aligned>);

template <uint64_t TSigma>
using Interleaved_16Aligned = GenericOccTable<rankvector::InterleavedBitvector16Aligned<TSigma>, "interleaved 16", "i16">;
static_assert(checkOccTable<Interleaved_16Aligned>);

template <uint64_t TSigma>
using Interleaved_32Aligned = GenericOccTable<rankvector::InterleavedBitvector32Aligned<TSigma>, "interleaved 32", "i32">;
static_assert(checkOccTable<Interleaved_32Aligned>);


template <uint64_t TSigma>
using Epr_8 = GenericOccTable<rankvector::InterleavedEPR8<TSigma>, "Epr (8bit)", "epr_8">;
static_assert(checkOccTable<Epr_8>);

template <uint64_t TSigma>
using Epr_16 = GenericOccTable<rankvector::InterleavedEPR16<TSigma>, "Epr (16bit)", "epr_16">;
static_assert(checkOccTable<Epr_16>);

template <uint64_t TSigma>
using Epr_32 = GenericOccTable<rankvector::InterleavedEPR32<TSigma>, "Epr (32bit)", "epr_32">;
static_assert(checkOccTable<Epr_32>);

template <uint64_t TSigma>
using Epr_8Aligned = GenericOccTable<rankvector::InterleavedEPR8Aligned<TSigma>, "Epr (8bit, aligned)", "epr_8a">;
static_assert(checkOccTable<Epr_8Aligned>);

template <uint64_t TSigma>
using Epr_16Aligned = GenericOccTable<rankvector::InterleavedEPR16Aligned<TSigma>, "Epr (16bit, aligned)", "epr_16a">;
static_assert(checkOccTable<Epr_16Aligned>);

template <uint64_t TSigma>
using Epr_32Aligned = GenericOccTable<rankvector::InterleavedEPR32Aligned<TSigma>, "Epr (32bit, aligned)", "epr_32a">;
static_assert(checkOccTable<Epr_32Aligned>);


template <uint64_t TSigma>
using EprV2_8 = GenericOccTable<rankvector::InterleavedEPRV2_8<TSigma>, "EprV2 (8bit)", "eprv2_8">;
static_assert(checkOccTable<EprV2_8>);

template <uint64_t TSigma>
using EprV2_16 = GenericOccTable<rankvector::InterleavedEPRV2_16<TSigma>, "EprV2 (16bit)", "eprv2_16">;
static_assert(checkOccTable<EprV2_16>);

template <uint64_t TSigma>
using EprV2_32 = GenericOccTable<rankvector::InterleavedEPRV2_32<TSigma>, "EprV2 (32bit)", "eprv2_32">;
static_assert(checkOccTable<EprV2_32>);

template <uint64_t TSigma>
using EprV2_8Aligned = GenericOccTable<rankvector::InterleavedEPRV2_8Aligned<TSigma>, "EprV2 (8bit, aligned)", "eprv2_8a">;
static_assert(checkOccTable<EprV2_8Aligned>);

template <uint64_t TSigma>
using EprV2_16Aligned = GenericOccTable<rankvector::InterleavedEPRV2_16Aligned<TSigma>, "EprV2 (16bit, aligned)", "eprv2_16a">;
static_assert(checkOccTable<EprV2_16Aligned>);

template <uint64_t TSigma>
using EprV2_32Aligned = GenericOccTable<rankvector::InterleavedEPRV2_32Aligned<TSigma>, "EprV2 (32bit, aligned)", "eprv2_32a">;
static_assert(checkOccTable<EprV2_32Aligned>);


template <uint64_t TSigma>
using EprV3_8 = GenericOccTable<rankvector::EPRV3_8<TSigma>, "EPRV3 (8bit)", "eprv3_8">;
static_assert(checkOccTable<EprV3_8>);

template <uint64_t TSigma>
using EprV3_16 = GenericOccTable<rankvector::EPRV3_16<TSigma>, "EPRV3 (16bit)", "eprv3_16">;
static_assert(checkOccTable<EprV3_16>);

template <uint64_t TSigma>
using EprV3_32 = GenericOccTable<rankvector::EPRV3_32<TSigma>, "EPRV3 (32bit)", "eprv3_32">;
static_assert(checkOccTable<EprV3_32>);


template <uint64_t TSigma>
using EprV4 = GenericOccTable<rankvector::EPRV4<TSigma>, "EPRV4", "eprv4">;
static_assert(checkOccTable<EprV4>);

template <uint64_t TSigma>
using EprV5 = GenericOccTable<rankvector::EPRV5<TSigma>, "EPRV5", "eprv5">;
static_assert(checkOccTable<EprV5>);

template <uint64_t TSigma>
using EprV6 = GenericOccTable<rankvector::DenseEPRV6<TSigma>, "EPRV6", "eprv6">;
static_assert(checkOccTable<EprV6>);

template <uint64_t TSigma>
using EprV7 = GenericOccTable<rankvector::InterleavedEPRV7<TSigma>, "EPRV7", "eprv7">;
static_assert(checkOccTable<EprV7>);

template <uint64_t TSigma>
using InterleavedWavelet = GenericOccTable<rankvector::InterleavedWavelet<TSigma>, "Interleaved Wavelet", "iw">;
static_assert(checkOccTable<InterleavedWavelet>);

template <uint64_t TSigma>
using Wavelet = GenericOccTable<rankvector::Wavelet<TSigma>, "Wavelet", "w">;
static_assert(checkOccTable<Wavelet>);

template <uint64_t TSigma>
using Sdsl_wt_bldc = GenericOccTable<rankvector::Sdsl_wt_bldc<TSigma>, "SDSL wt_bldc", "sdsl_wt_bldc">;
static_assert(checkOccTable<Sdsl_wt_bldc>);

}
