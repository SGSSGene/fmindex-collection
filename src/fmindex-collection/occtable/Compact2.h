#pragma once

#include "GenericInterleaved.h"

namespace fmindex_collection {
namespace occtable {
namespace compact2 {

template <size_t TSigma>
struct OccTable : genericInterleaved::OccTable<TSigma, 8, uint16_t> {
    static auto name() -> std::string {
        return "Interleaved 16bit";
    }

    static auto extension() -> std::string {
        return "i16";
    }
};

}
}
}
