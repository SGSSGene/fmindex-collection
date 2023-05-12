// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file.
// -----------------------------------------------------------------------------------------------------
#pragma once

#include "GenericInterleaved.h"

namespace fmindex_collection {
namespace occtable {

namespace interleaved8 {
template <uint64_t TSigma>
struct OccTable : genericInterleaved::OccTable<TSigma, 8, uint8_t> {
    using genericInterleaved::OccTable<TSigma, 8, uint8_t>::OccTable;
    static auto name() -> std::string {
        return "Interleaved 8bit";
    }

    static auto extension() -> std::string {
        return "i8";
    }
};
static_assert(checkOccTable<OccTable>);
}

namespace interleaved16 {
template <uint64_t TSigma>
struct OccTable : genericInterleaved::OccTable<TSigma, 8, uint16_t> {
    using genericInterleaved::OccTable<TSigma, 8, uint16_t>::OccTable;
    static auto name() -> std::string {
        return "Interleaved 16bit";
    }

    static auto extension() -> std::string {
        return "i16";
    }
};
static_assert(checkOccTable<OccTable>);
}

namespace interleaved32 {
template <uint64_t TSigma>
struct OccTable : genericInterleaved::OccTable<TSigma, 8, uint32_t> {
    using genericInterleaved::OccTable<TSigma, 8, uint32_t>::OccTable;
    static auto name() -> std::string {
        return "Interleaved 32bit";
    }

    static auto extension() -> std::string {
        return "i32";
    }
};
static_assert(checkOccTable<OccTable>);
}

namespace interleaved8Aligned {
template <uint64_t TSigma>
struct OccTable : genericInterleaved::OccTable<TSigma, 64, uint8_t> {
    using genericInterleaved::OccTable<TSigma, 64, uint8_t>::OccTable;
    static auto name() -> std::string {
        return "Aligned Interleaved 8bit";
    }

    static auto extension() -> std::string {
        return "i8a";
    }
};
static_assert(checkOccTable<OccTable>);
}

namespace interleaved16Aligned {
template <uint64_t TSigma>
struct OccTable : genericInterleaved::OccTable<TSigma, 64, uint16_t> {
    using genericInterleaved::OccTable<TSigma, 64, uint16_t>::OccTable;
    static auto name() -> std::string {
        return "Aligned Interleaved 16bit";
    }

    static auto extension() -> std::string {
        return "i16a";
    }
};
static_assert(checkOccTable<OccTable>);
}

namespace interleaved32Aligned {
template <uint64_t TSigma>
struct OccTable : genericInterleaved::OccTable<TSigma, 64, uint32_t> {
    using genericInterleaved::OccTable<TSigma, 64, uint32_t>::OccTable;

    static auto name() -> std::string {
        return "Aligned Interleaved 32bit";
    }

    static auto extension() -> std::string {
        return "i32a";
    }
};
static_assert(checkOccTable<OccTable>);
}



}
}
