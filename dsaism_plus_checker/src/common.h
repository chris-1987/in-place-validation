#ifndef __COMMON_H
#define __COMMON_H

#include "stxxl/bits/common/uint_types.h"
#include "stxxl/vector"
#include "stxxl/bits/containers/queue.h"
#include "stxxl/bits/containers/deque.h"
#include "stxxl/bits/containers/priority_queue.h"
#include "stxxl/bits/containers/sorter.h"
#include "stxxl/bits/io/syscall_file.h"

#include <limits>
#include <vector>
#include <utility>

// alias for integral types
using uint8 = stxxl::uint8;

using uint16 = stxxl::uint16;

using uint32 = stxxl::uint32;

using uint40 = stxxl::uint40;

using uint64 = stxxl::uint64;

// L-type or S-type
constexpr uint8 L_TYPE = 0;

constexpr uint8 S_TYPE = 1;

// LONG_TYPE or SHORT_TYPE
constexpr uint8 SHORT_TYPE = 0;

constexpr uint8 LONG_TYPE = 1;

// D_LOW & D_HIGH
constexpr uint8 D_LOW = 4;

constexpr uint8 D_HIGH = 4;

// memory allocation
constexpr uint64 K_512 = 512 * 1024;

constexpr uint64 K_1024 = 1024 * 1024;

constexpr uint64 MAX_MEM = 3 * 1024 * K_1024;

constexpr uint64 MAX_MEM2 = MAX_MEM;

constexpr uint64 MAX_ITEM = 16 * K_1024;




#endif // common_h
