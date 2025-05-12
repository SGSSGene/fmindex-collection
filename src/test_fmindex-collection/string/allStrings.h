// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#pragma once

#include <fmindex-collection/bitvector/all.h>
#include <fmindex-collection/string/all.h>

#define ALLRANKVECTORS_IMPL(Sigma) \
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
    fmindex_collection::string::PairedL0L1_NEPRV9_2048_64k<Sigma>

#if FMC_USE_SDSL
#define ALLRANKVECTORS(Sigma) \
    ALLRANKVECTORS_IMPL(Sigma), \
    fmindex_collection::string::Sdsl_wt_bldc<Sigma>, \
    fmindex_collection::string::Sdsl_wt_epr<Sigma>
#else
#define ALLRANKVECTORS(Sigma) ALLRANKVECTORS_IMPL(Sigma)
#endif

#define ALLLARGERANKVECTORS(Sigma) \
    fmindex_collection::string::L0L1_NEPRV9_64_64k<Sigma>, \
    fmindex_collection::string::L0L1_NEPRV9_512_64k<Sigma>, \
    fmindex_collection::string::L0L1_NEPRV9_4096_64k<Sigma>, \
    fmindex_collection::string::PairedL0L1_NEPRV9_64_64k<Sigma>, \
    fmindex_collection::string::PairedL0L1_NEPRV9_512_64k<Sigma>, \
    fmindex_collection::string::PairedL0L1_NEPRV9_4096_64k<Sigma>
