#pragma once

#include "GenericInterleaved.h"

namespace occtable {
namespace compactAligned {

template <size_t TSigma>
struct OccTable : genericInterleaved::OccTable<TSigma, 64, uint32_t> {
    static auto name() -> std::string {
        return "Interleaved 32bit Aligned";
    }

    static auto extension() -> std::string {
        return "i32a";
    }
};

}
}
