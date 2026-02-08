// SPDX-FileCopyrightText: 2026 Simon Gene Gottlieb
// SPDX-License-Identifier: CC0-1.0
#include "allStrings.h"
#include "utils.h"

#include <fmindex-collection/string/PairedFlattenedBitvectors2L_b.h>
#include <fmindex-collection/string/AdapterStringKStep.h>


TEST_CASE("check if rank on strings with 'dual_limit' functions work", "[string][string_limit]") {
    using T1 = fmc::string::PairedFlattenedBitvectors2L_b<16, 2, 512, 65536>;
    static_assert(fmc::StringKStep_c<T1>);

    using T2 = fmc::string::InterleavedBitvector16<16>;
    static_assert(!fmc::StringKStep_c<T2>);

    using T3 = fmc::string::AdapterStringKStep<16, 2, fmc::string::InterleavedBitvector16>;
    static_assert(fmc::StringKStep_c<T3>);
}
