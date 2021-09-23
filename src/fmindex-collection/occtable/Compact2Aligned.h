#pragma once

#include "GenericInterleaved.h"

namespace occtable {
namespace compact2Aligned {

template <size_t TSigma>
struct OccTable : genericInterleaved::OccTable<TSigma, 64, uint16_t> {};

}
}
