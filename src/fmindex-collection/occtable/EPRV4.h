// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "GenericOccTable.h"
#include "../rankvector/EPRV4.h"

namespace fmindex_collection {
namespace occtable {
template <uint64_t TSigma>
using EprV4 = GenericOccTable<rankvector::EPRV4<TSigma>, "EPRV4", "eprv4">;

static_assert(checkOccTable<EprV4>);

}
}
