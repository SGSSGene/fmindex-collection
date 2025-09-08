// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#pragma once

#include <fmindex-collection/bitvector/all.h>
#include <fmindex-collection/bitvector/RBBitvector.h>

#define ALLBITVECTORS \
    fmc::bitvector::L0_64Bitvector, \
    fmc::bitvector::L0_128Bitvector, \
    fmc::bitvector::L0_256Bitvector, \
    fmc::bitvector::L0_512Bitvector, \
    fmc::bitvector::L0_1024Bitvector, \
    fmc::bitvector::L0_2048Bitvector, \
    fmc::bitvector::PairedL0_64Bitvector, \
    fmc::bitvector::PairedL0_128Bitvector, \
    fmc::bitvector::PairedL0_256Bitvector, \
    fmc::bitvector::PairedL0_512Bitvector, \
    fmc::bitvector::PairedL0_1024Bitvector, \
    fmc::bitvector::PairedL0_2048Bitvector, \
    fmc::bitvector::L0L1_64_64kBitvector, \
    fmc::bitvector::L0L1_128_64kBitvector, \
    fmc::bitvector::L0L1_256_64kBitvector, \
    fmc::bitvector::L0L1_512_64kBitvector, \
    fmc::bitvector::L0L1_1024_64kBitvector, \
    fmc::bitvector::L0L1_2048_64kBitvector, \
    fmc::bitvector::PairedL0L1_64_64kBitvector, \
    fmc::bitvector::PairedL0L1_128_64kBitvector, \
    fmc::bitvector::PairedL0L1_256_64kBitvector, \
    fmc::bitvector::PairedL0L1_512_64kBitvector, \
    fmc::bitvector::PairedL0L1_1024_64kBitvector, \
    fmc::bitvector::PairedL0L1_2048_64kBitvector

#define ALLSPARSEBITVECTORS \
    fmc::bitvector::L0_64Bitvector, \
    fmc::bitvector::L0L1_512_64kBitvector, \
    fmc::bitvector::RBBitvector<2, fmc::bitvector::L0L1_512_64kBitvector, fmc::bitvector::L0L1_512_64kBitvector>, \
    fmc::bitvector::RBBitvector<2, fmc::bitvector::L0L1_512_64kBitvector, fmc::bitvector::L0_64Bitvector>, \
    fmc::bitvector::SparseBLEBitvector<2, fmc::bitvector::L0_64Bitvector, fmc::bitvector::L0_64Bitvector>, \
    fmc::bitvector::SparseBLEBitvector<2, fmc::bitvector::L0_64Bitvector, fmc::bitvector::L0L1_512_64kBitvector>
