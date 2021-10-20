#pragma once

#include "GenericInterleaved.h"

namespace fmindex_collection {
namespace occtable {
namespace compact {

template <size_t TSigma>
struct OccTable : genericInterleaved::OccTable<TSigma, 8, uint32_t> {
    static auto name() -> std::string {
        return "Interleaved 32bit";
    }

    static auto extension() -> std::string {
        return "i32";
    }
};

}
}
}
