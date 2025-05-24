// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#pragma once

#include <fmindex-collection/bitvector/all.h>
#include <fmindex-collection/string/all.h>

#define ALLSTRINGSWITHRANK_IMPL(Sigma) \
    fmindex_collection::string::InterleavedBitvector16<Sigma>, \
    fmindex_collection::string::L0L1_NEPRV9_64_64k<Sigma>, \
    fmindex_collection::string::L0L1_NEPRV9_128_64k<Sigma>, \
    fmindex_collection::string::L0L1_NEPRV9_256_64k<Sigma>, \
    fmindex_collection::string::L0L1_NEPRV9_512_64k<Sigma>, \
    fmindex_collection::string::L0L1_NEPRV9_1024_64k<Sigma>, \
    fmindex_collection::string::L0L1_NEPRV9_2048_64k<Sigma>, \
    fmindex_collection::string::PairedL0L1_NEPRV9_64_64k<Sigma>, \
    fmindex_collection::string::PairedL0L1_NEPRV9_128_64k<Sigma>, \
    fmindex_collection::string::PairedL0L1_NEPRV9_256_64k<Sigma>, \
    fmindex_collection::string::PairedL0L1_NEPRV9_512_64k<Sigma>, \
    fmindex_collection::string::PairedL0L1_NEPRV9_1024_64k<Sigma>, \
    fmindex_collection::string::PairedL0L1_NEPRV9_2048_64k<Sigma>, \
    fmindex_collection::string::MultiBitvector_Bitvector<Sigma>, \
    fmindex_collection::string::InterleavedEPR16<Sigma>, \
    fmindex_collection::string::InterleavedEPRV2_16<Sigma>, \
    fmindex_collection::string::MultiaryWavelet<Sigma, fmindex_collection::string::PairedL0L1_NEPRV9_64_64k>, \
    fmindex_collection::string::MultiaryWavelet<Sigma>, \
    fmindex_collection::string::MultiaryWavelet<Sigma, fmindex_collection::string::PairedL0L1_NEPRV9_512_64k, 4, fmindex_collection::string::MultiaryWavelet>, \
    fmindex_collection::string::MultiaryWavelet<Sigma, fmindex_collection::string::PairedL0L1_NEPRV9_512_64k, 8, fmindex_collection::string::MultiaryWavelet>


#if FMC_USE_SDSL
#define ALLSTRINGSWITHRANK(Sigma) \
    ALLSTRINGSWITHRANK_IMPL(Sigma), \
    fmindex_collection::string::Sdsl_wt_bldc<Sigma>, \
    fmindex_collection::string::Sdsl_wt_epr<Sigma>
#else
#define ALLSTRINGSWITHRANK(Sigma) ALLSTRINGSWITHRANK_IMPL(Sigma)
#endif

#define ALLLARGESTRINGSWITHRANK(Sigma) \
    /*fmindex_collection::string::L0L1_NEPRV9_512_64k<Sigma>,*/ \
    /*fmindex_collection::string::L0L1_NEPRV9_2048_64k<Sigma>,*/ \
    /*fmindex_collection::string::PairedL0L1_NEPRV9_512_64k<Sigma>,*/ \
    /*fmindex_collection::string::PairedL0L1_NEPRV9_2048_64k<Sigma>,*/ \
    fmindex_collection::string::MultiaryWavelet<Sigma, fmindex_collection::string::PairedL0L1_NEPRV9_64_64k>, \
    fmindex_collection::string::MultiaryWavelet<Sigma>, \
    fmindex_collection::string::MultiaryWavelet<Sigma, fmindex_collection::string::PairedL0L1_NEPRV9_512_64k, 4, fmindex_collection::string::MultiaryWavelet>, \
    fmindex_collection::string::MultiaryWavelet<Sigma, fmindex_collection::string::PairedL0L1_NEPRV9_512_64k, 8, fmindex_collection::string::MultiaryWavelet>, \
    fmindex_collection::string::MultiaryWavelet<Sigma, fmindex_collection::string::PairedL0L1_NEPRV9_512_64k, 16, fmindex_collection::string::MultiaryWavelet>, \
    fmindex_collection::string::MultiaryWavelet<Sigma, fmindex_collection::string::PairedL0L1_NEPRV9_512_64k, 32, fmindex_collection::string::MultiaryWavelet>, \
    fmindex_collection::string::MultiaryWavelet<Sigma, fmindex_collection::string::PairedL0L1_NEPRV9_512_64k, 64, fmindex_collection::string::MultiaryWavelet>
