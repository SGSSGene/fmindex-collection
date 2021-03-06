#pragma once

#include "Bitvector.h"
#include "BitvectorPrefix.h"
#include "InterleavedBitvector.h"
#include "InterleavedPrefix.h"
#include "InterleavedWavelet.h"
#include "InterleavedWavelet32.h"
#include "InterleavedEPR.h"
#include "InterleavedEPRV2.h"
#include "EPRV3.h"
#include "EPRV4.h"
#include "EPRV5.h"
#include "DenseEPRV6.h"
#include "Naive.h"
#include "Wavelet.h"

#ifdef FMC_USE_SDSL
#include "Sdsl_wt_bldc.h"
#include "Sdsl_wt_epr.h"
#endif
