// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0
#pragma once

#include <fmindex-collection/bitvector/all.h>
#include <fmindex-collection/string/all.h>

#define ALLSTRINGSWITHRANK_IMPL \
    fmc::string::InterleavedBitvector16, \
    fmc::string::FlattenedBitvectors_64_64k, \
    fmc::string::FlattenedBitvectors_128_64k, \
    fmc::string::FlattenedBitvectors_256_64k, \
    fmc::string::FlattenedBitvectors_512_64k, \
    fmc::string::FlattenedBitvectors_1024_64k, \
    fmc::string::FlattenedBitvectors_2048_64k, \
    fmc::string::PairedFlattenedBitvectors_64_64k, \
    fmc::string::PairedFlattenedBitvectors_128_64k, \
    fmc::string::PairedFlattenedBitvectors_256_64k, \
    fmc::string::PairedFlattenedBitvectors_512_64k, \
    fmc::string::PairedFlattenedBitvectors_1024_64k, \
    fmc::string::PairedFlattenedBitvectors_2048_64k, \
    fmc::string::MultiBitvector_Bitvector, \
    fmc::string::InterleavedEPR16, \
    fmc::string::InterleavedEPRV2_16, \
    fmc::string::MultiaryWavelet_64_64k, \
    fmc::string::MultiaryWavelet_512_64k, \
    fmc::string::MultiaryWavelet_s16, \
    fmc::string::MultiaryWavelet_s256


#if FMC_USE_SDSL
#define ALLSTRINGSWITHRANK \
    ALLSTRINGSWITHRANK_IMPL, \
    fmc::string::Sdsl_wt_bldc, \
    fmc::string::Sdsl_wt_epr
#else
#define ALLSTRINGSWITHRANK ALLSTRINGSWITHRANK_IMPL
#endif
#if 1
#define ALLLARGESTRINGSWITHRANK \
    /*fmc::string::FlattenedBitvectors_512_64k,*/ \
    /*fmc::string::FlattenedBitvectors_2048_64k,*/ \
    /*fmc::string::PairedFlattenedBitvectors_512_64k,*/ \
    /*fmc::string::PairedFlattenedBitvectors_2048_64k,*/ \
    fmc::string::MultiaryWavelet<0, fmc::string::PairedFlattenedBitvectors_64_64k>::partial_t, \
    fmc::string::MultiaryWavelet<0>::partial_t, \
    fmc::string::MultiaryWavelet<0, fmc::string::PairedFlattenedBitvectors_512_64k, fmc::string::MultiaryWavelet, 4>::partial_t, \
    fmc::string::MultiaryWavelet<0, fmc::string::PairedFlattenedBitvectors_512_64k, fmc::string::MultiaryWavelet, 8>::partial_t, \
    fmc::string::MultiaryWavelet<0, fmc::string::PairedFlattenedBitvectors_512_64k, fmc::string::MultiaryWavelet, 16>::partial_t, \
    fmc::string::MultiaryWavelet<0, fmc::string::PairedFlattenedBitvectors_512_64k, fmc::string::MultiaryWavelet, 32>::partial_t, \
    fmc::string::MultiaryWavelet<0, fmc::string::PairedFlattenedBitvectors_512_64k, fmc::string::MultiaryWavelet, 64>::partial_t
#else
#define ALLLARGESTRINGSWITHRANK
#endif

