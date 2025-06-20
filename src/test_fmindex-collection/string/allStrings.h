// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#pragma once

#include <fmindex-collection/bitvector/all.h>
#include <fmindex-collection/string/all.h>

#define ALLSTRINGSWITHRANK_IMPL \
    fmindex_collection::string::InterleavedBitvector16, \
    fmindex_collection::string::FlattenedBitvectors_64_64k, \
    fmindex_collection::string::FlattenedBitvectors_128_64k, \
    fmindex_collection::string::FlattenedBitvectors_256_64k, \
    fmindex_collection::string::FlattenedBitvectors_512_64k, \
    fmindex_collection::string::FlattenedBitvectors_1024_64k, \
    fmindex_collection::string::FlattenedBitvectors_2048_64k, \
    fmindex_collection::string::PairedFlattenedBitvectors_64_64k, \
    fmindex_collection::string::PairedFlattenedBitvectors_128_64k, \
    fmindex_collection::string::PairedFlattenedBitvectors_256_64k, \
    fmindex_collection::string::PairedFlattenedBitvectors_512_64k, \
    fmindex_collection::string::PairedFlattenedBitvectors_1024_64k, \
    fmindex_collection::string::PairedFlattenedBitvectors_2048_64k, \
    fmindex_collection::string::MultiBitvector_Bitvector, \
    fmindex_collection::string::InterleavedEPR16, \
    fmindex_collection::string::InterleavedEPRV2_16, \
    fmindex_collection::string::MultiaryWavelet_64_64k, \
    fmindex_collection::string::MultiaryWavelet_512_64k, \
    fmindex_collection::string::MultiaryWavelet_s16, \
    fmindex_collection::string::MultiaryWavelet_s256


#if FMC_USE_SDSL
#define ALLSTRINGSWITHRANK \
    ALLSTRINGSWITHRANK_IMPL, \
    fmindex_collection::string::Sdsl_wt_bldc, \
    fmindex_collection::string::Sdsl_wt_epr
#else
#define ALLSTRINGSWITHRANK ALLSTRINGSWITHRANK_IMPL
#endif
#if 1
#define ALLLARGESTRINGSWITHRANK \
    /*fmindex_collection::string::FlattenedBitvectors_512_64k,*/ \
    /*fmindex_collection::string::FlattenedBitvectors_2048_64k,*/ \
    /*fmindex_collection::string::PairedFlattenedBitvectors_512_64k,*/ \
    /*fmindex_collection::string::PairedFlattenedBitvectors_2048_64k,*/ \
    fmindex_collection::string::MultiaryWavelet<0, fmindex_collection::string::PairedFlattenedBitvectors_64_64k>::partial_t, \
    fmindex_collection::string::MultiaryWavelet<0>::partial_t, \
    fmindex_collection::string::MultiaryWavelet<0, fmindex_collection::string::PairedFlattenedBitvectors_512_64k, fmindex_collection::string::MultiaryWavelet, 4>::partial_t, \
    fmindex_collection::string::MultiaryWavelet<0, fmindex_collection::string::PairedFlattenedBitvectors_512_64k, fmindex_collection::string::MultiaryWavelet, 8>::partial_t, \
    fmindex_collection::string::MultiaryWavelet<0, fmindex_collection::string::PairedFlattenedBitvectors_512_64k, fmindex_collection::string::MultiaryWavelet, 16>::partial_t, \
    fmindex_collection::string::MultiaryWavelet<0, fmindex_collection::string::PairedFlattenedBitvectors_512_64k, fmindex_collection::string::MultiaryWavelet, 32>::partial_t, \
    fmindex_collection::string::MultiaryWavelet<0, fmindex_collection::string::PairedFlattenedBitvectors_512_64k, fmindex_collection::string::MultiaryWavelet, 64>::partial_t
#else
#define ALLLARGESTRINGSWITHRANK
#endif

