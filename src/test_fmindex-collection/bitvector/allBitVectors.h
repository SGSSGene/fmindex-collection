// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#pragma once

#include <fmindex-collection/bitvector/all.h>


#if !__clang__ || __clang_major__ >= 19 // !TODO workaround (weird optimization bug?)
#define AddIfNotClang18OrOlder(x) x,
#else
#define AddIfNotClang18OrOlder(x)
#endif

#define ALLBITVECTORS \
    fmindex_collection::bitvector::Bitvector, \
    fmindex_collection::bitvector::PrunedBitvector, \
    fmindex_collection::bitvector::CompactBitvector, \
    fmindex_collection::bitvector::CompactBitvector4Blocks, \
    fmindex_collection::bitvector::L1Bitvector, \
    fmindex_collection::bitvector::L1_64Bitvector, \
    fmindex_collection::bitvector::L1_128Bitvector, \
    fmindex_collection::bitvector::L1_256Bitvector, \
    fmindex_collection::bitvector::L1_512Bitvector, \
    fmindex_collection::bitvector::DoubleL1_64Bitvector, \
    fmindex_collection::bitvector::DoubleL1_128Bitvector, \
    fmindex_collection::bitvector::DoubleL1_256Bitvector, \
    fmindex_collection::bitvector::DoubleL1_512Bitvector, \
    fmindex_collection::bitvector::DoubleL1L2_64_4kBitvector, \
    fmindex_collection::bitvector::DoubleL1L2_128_4kBitvector, \
    fmindex_collection::bitvector::DoubleL1L2_256_4kBitvector, \
    fmindex_collection::bitvector::DoubleL1L2_512_4kBitvector, \
    fmindex_collection::bitvector::DoubleL1L2_64_64kBitvector, \
    fmindex_collection::bitvector::DoubleL1L2_128_64kBitvector, \
    fmindex_collection::bitvector::DoubleL1L2_256_64kBitvector, \
    fmindex_collection::bitvector::DoubleL1L2_2048_64kBitvector, \
    AddIfNotClang18OrOlder(fmindex_collection::bitvector::CompactDoubleL1L2_64_4kBitvector) \
    AddIfNotClang18OrOlder(fmindex_collection::bitvector::CompactDoubleL1L2_128_4kBitvector) \
    AddIfNotClang18OrOlder(fmindex_collection::bitvector::CompactDoubleL1L2_256_4kBitvector) \
    AddIfNotClang18OrOlder(fmindex_collection::bitvector::CompactDoubleL1L2_512_4kBitvector) \
    AddIfNotClang18OrOlder(fmindex_collection::bitvector::CompactDoubleL1L2_64_64kBitvector) \
    AddIfNotClang18OrOlder(fmindex_collection::bitvector::CompactDoubleL1L2_128_64kBitvector) \
    AddIfNotClang18OrOlder(fmindex_collection::bitvector::CompactDoubleL1L2_256_64kBitvector) \
    AddIfNotClang18OrOlder(fmindex_collection::bitvector::CompactDoubleL1L2_2048_64kBitvector) \
    fmindex_collection::bitvector::SparseBLEBitvector<2>, \
    fmindex_collection::bitvector::SparseBLEBitvector<-2>, \
    (fmindex_collection::bitvector::SparseBLEBitvector<3, fmindex_collection::bitvector::SparseBLEBitvector<2>>), \
    (fmindex_collection::bitvector::SparseBLEBitvector<3, fmindex_collection::bitvector::Bitvector, fmindex_collection::bitvector::SparseBLEBitvector<2>>), \
    (fmindex_collection::bitvector::SparseBLEBitvector<2, fmindex_collection::bitvector::DoubleL1L2_64_4kBitvector, fmindex_collection::bitvector::DoubleL1L2_64_4kBitvector>)


//(fmindex_collection::bitvector::DoubleL1L2_NBitvector<4, 16>),
