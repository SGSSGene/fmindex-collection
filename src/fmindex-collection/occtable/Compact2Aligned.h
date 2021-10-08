#pragma once

#include "GenericInterleaved.h"

namespace occtable {
namespace compact2Aligned {

template <size_t TSigma>
struct OccTable : genericInterleaved::OccTable<TSigma, 64, uint16_t> {
    static auto name() -> std::string {
        return "Interleaved 16bit Aligned";
    }

    static auto extension() -> std::string {
        return "i16a";
    }
};

}
}
