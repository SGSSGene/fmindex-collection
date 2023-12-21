// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause
#pragma once

#include "GenericOccTable.h"
#include "../rankvector/Naive.h"

namespace fmindex_collection {
namespace occtable {
template <uint64_t TSigma>
using Naive = GenericOccTable<rankvector::Naive<TSigma>, "Naive", "n">;

static_assert(checkOccTable<Naive>);

}
}
