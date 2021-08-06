#pragma once

#include "GenericInterleaved.h"

namespace occtable {
namespace compact {

template <size_t TSigma>
struct OccTable : genericInterleaved::OccTable<TSigma, 8, uint32_t> {};

}
}
