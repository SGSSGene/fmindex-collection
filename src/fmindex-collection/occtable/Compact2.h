#pragma once

#include "GenericInterleaved.h"

namespace occtable {
namespace compact2 {

template <size_t TSigma>
struct OccTable : genericInterleaved::OccTable<TSigma, 8, uint16_t> {};

}
}
