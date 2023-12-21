// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "GenericOccTable.h"
#include "../rankvector/EPRV3.h"

namespace fmindex_collection {
namespace occtable {
template <uint64_t TSigma>
using EprV3_8 = GenericOccTable<rankvector::EPRV3_8<TSigma>, "EPRV3 (8bit)", "eprv3_8">;

template <uint64_t TSigma>
using EprV3_16 = GenericOccTable<rankvector::EPRV3_16<TSigma>, "EPRV3 (16bit)", "eprv3_16">;

template <uint64_t TSigma>
using EprV3_32 = GenericOccTable<rankvector::EPRV3_32<TSigma>, "EPRV3 (32bit)", "eprv3_32">;

static_assert(checkOccTable<EprV3_8>);
static_assert(checkOccTable<EprV3_16>);
static_assert(checkOccTable<EprV3_32>);

}
}
