// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#pragma once

#include <fmindex-collection/bitvector/all.h>
#define ALLBITVECTORS \
    fmindex_collection::bitvector::Bitvector, \
    fmindex_collection::bitvector::CompactBitvector, \
    fmindex_collection::bitvector::CompactBitvector4Blocks, \
    fmindex_collection::bitvector::L1Bitvector, \
    fmindex_collection::bitvector::L1_64Bitvector, \
    fmindex_collection::bitvector::L1_128Bitvector, \
    fmindex_collection::bitvector::L1_256Bitvector, \
    fmindex_collection::bitvector::DoubleL1_64Bitvector, \
    fmindex_collection::bitvector::DoubleL1_128Bitvector, \
    fmindex_collection::bitvector::DoubleL1_256Bitvector, \
    fmindex_collection::bitvector::SparseBLEBitvector<2>, \
    fmindex_collection::bitvector::SparseBLEBitvector<-2>, \
    (fmindex_collection::bitvector::SparseBLEBitvector<3, fmindex_collection::bitvector::SparseBLEBitvector<2>>), \
    (fmindex_collection::bitvector::SparseBLEBitvector<3, fmindex_collection::bitvector::Bitvector, fmindex_collection::bitvector::SparseBLEBitvector<2>>)
