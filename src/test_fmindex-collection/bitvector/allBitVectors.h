// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#pragma once

#include <fmindex-collection/bitvector/all.h>

#ifdef FMC_USE_PASTA
#include "Pasta_FlatRank.h"
#include "Pasta_WideRank.h"
#endif

#ifdef FMC_USE_SDSL
#include "sdsl_v.h"
#include "sdsl_v5.h"
#endif

#ifdef FMC_USE_SUX
#include "sux_Rank9.h"
#endif

#ifdef FMC_USE_RANKSELECT
#include "RankSelect_Rank.h"
#endif



#ifndef __EMSCRIPTEN__
#define ALLBITVECTORS \
    fmc::bitvector::Bitvector1L_64, \
    fmc::bitvector::Bitvector1L_128, \
    fmc::bitvector::Bitvector1L_256, \
    fmc::bitvector::Bitvector1L_512, \
    fmc::bitvector::Bitvector1L_1024, \
    fmc::bitvector::Bitvector1L_2048, \
    fmc::bitvector::PairedBitvector1L_64, \
    fmc::bitvector::PairedBitvector1L_128, \
    fmc::bitvector::PairedBitvector1L_256, \
    fmc::bitvector::PairedBitvector1L_512, \
    fmc::bitvector::PairedBitvector1L_1024, \
    fmc::bitvector::PairedBitvector1L_2048, \
    fmc::bitvector::Bitvector2L_64_64k, \
    fmc::bitvector::Bitvector2L_128_64k, \
    fmc::bitvector::Bitvector2L_256_64k, \
    fmc::bitvector::Bitvector2L_512_64k, \
    fmc::bitvector::Bitvector2L_1024_64k, \
    fmc::bitvector::Bitvector2L_2048_64k, \
    fmc::bitvector::PairedBitvector2L_64_64k, \
    fmc::bitvector::PairedBitvector2L_128_64k, \
    fmc::bitvector::PairedBitvector2L_256_64k, \
    fmc::bitvector::PairedBitvector2L_512_64k, \
    fmc::bitvector::PairedBitvector2L_1024_64k, \
    fmc::bitvector::PairedBitvector2L_2048_64k
#else
#define ALLBITVECTORS \
    fmc::bitvector::Bitvector1L_64, \
    fmc::bitvector::Bitvector1L_512, \
    fmc::bitvector::PairedBitvector1L_64, \
    fmc::bitvector::PairedBitvector1L_512, \
    fmc::bitvector::Bitvector2L_64_64k, \
    fmc::bitvector::Bitvector2L_512_64k, \
    fmc::bitvector::PairedBitvector2L_64_64k, \
    fmc::bitvector::PairedBitvector2L_512_64k
#endif

#define ALLSPARSEBITVECTORS \
    fmc::bitvector::Bitvector1L_64, \
    fmc::bitvector::Bitvector2L_512_64k, \
    fmc::bitvector::RBBitvector<2, fmc::bitvector::Bitvector2L_64_64k, fmc::bitvector::Bitvector1L_64>, \
    fmc::bitvector::RBBitvector<2, fmc::bitvector::Bitvector2L_512_64k, fmc::bitvector::Bitvector1L_64>, \
    fmc::bitvector::RBBitvector<2, fmc::bitvector::Bitvector2L_512_64k, fmc::bitvector::Bitvector2L_512_64k>, \
    fmc::bitvector::SparseBLEBitvector<2, fmc::bitvector::Bitvector1L_64, fmc::bitvector::Bitvector1L_64>, \
    fmc::bitvector::SparseBLEBitvector<2, fmc::bitvector::Bitvector1L_64, fmc::bitvector::Bitvector2L_512_64k>, \
    fmc::bitvector::SparseBLEBitvector<4, fmc::bitvector::Bitvector1L_64, fmc::bitvector::Bitvector1L_64>, \
    fmc::bitvector::SparseBLEBitvector<4, fmc::bitvector::Bitvector1L_64, fmc::bitvector::Bitvector2L_512_64k>, \
    fmc::bitvector::SparseDynRBBitvector<fmc::bitvector::Bitvector2L_64_64k, fmc::bitvector::Bitvector1L_64>, \
    fmc::bitvector::SparseDynRBBitvector<fmc::bitvector::Bitvector2L_512_64k, fmc::bitvector::Bitvector1L_64>, \
    fmc::bitvector::SparseDynRBBitvector<fmc::bitvector::Bitvector2L_512_64k, fmc::bitvector::Bitvector2L_512_64k>
