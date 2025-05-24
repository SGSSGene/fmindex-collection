// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#pragma once

#include <fmindex-collection/bitvector/all.h>

#define ALLBITVECTORS \
    fmindex_collection::bitvector::L0_64Bitvector, \
    fmindex_collection::bitvector::L0_128Bitvector, \
    fmindex_collection::bitvector::L0_256Bitvector, \
    fmindex_collection::bitvector::L0_512Bitvector, \
    fmindex_collection::bitvector::L0_1024Bitvector, \
    fmindex_collection::bitvector::L0_2048Bitvector, \
    fmindex_collection::bitvector::PairedL0_64Bitvector, \
    fmindex_collection::bitvector::PairedL0_128Bitvector, \
    fmindex_collection::bitvector::PairedL0_256Bitvector, \
    fmindex_collection::bitvector::PairedL0_512Bitvector, \
    fmindex_collection::bitvector::PairedL0_1024Bitvector, \
    fmindex_collection::bitvector::PairedL0_2048Bitvector, \
    fmindex_collection::bitvector::L0L1_64_64kBitvector, \
    fmindex_collection::bitvector::L0L1_128_64kBitvector, \
    fmindex_collection::bitvector::L0L1_256_64kBitvector, \
    fmindex_collection::bitvector::L0L1_512_64kBitvector, \
    fmindex_collection::bitvector::L0L1_1024_64kBitvector, \
    fmindex_collection::bitvector::L0L1_2048_64kBitvector, \
    fmindex_collection::bitvector::PairedL0L1_64_64kBitvector, \
    fmindex_collection::bitvector::PairedL0L1_128_64kBitvector, \
    fmindex_collection::bitvector::PairedL0L1_256_64kBitvector, \
    fmindex_collection::bitvector::PairedL0L1_512_64kBitvector, \
    fmindex_collection::bitvector::PairedL0L1_1024_64kBitvector, \
    fmindex_collection::bitvector::PairedL0L1_2048_64kBitvector

#define ALLSPARSEBITVECTORS \
    fmindex_collection::bitvector::PairedL0L1_512_64kBitvector, \
    fmindex_collection::bitvector::SparseBLEBitvector<1, fmindex_collection::bitvector::PairedL0L1_512_64kBitvector, fmindex_collection::bitvector::PairedL0L1_512_64kBitvector>, \
    fmindex_collection::bitvector::SparseBLEBitvector<2, fmindex_collection::bitvector::PairedL0L1_512_64kBitvector, fmindex_collection::bitvector::PairedL0L1_512_64kBitvector>, \
    fmindex_collection::bitvector::SparseBLEBitvector<3, fmindex_collection::bitvector::PairedL0L1_512_64kBitvector, fmindex_collection::bitvector::PairedL0L1_512_64kBitvector>, \
    fmindex_collection::bitvector::SparseBLEBitvector<4, fmindex_collection::bitvector::PairedL0L1_512_64kBitvector, fmindex_collection::bitvector::PairedL0L1_512_64kBitvector>, \
    fmindex_collection::bitvector::SparseBLEBitvector<5, fmindex_collection::bitvector::PairedL0L1_512_64kBitvector, fmindex_collection::bitvector::PairedL0L1_512_64kBitvector>
