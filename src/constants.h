#pragma once

#include <cstdint>
#include <algorithm>

#include "smallest_fit.h"

// Defines the extender itself. Each may be overridden from the build system,
// e.g. -DSNAPERZ_LENGTH=49, which is convenient for benchmarking and for the
// test harness that sweeps over many lengths in a single run.
#ifndef SNAPERZ_LENGTH
#define SNAPERZ_LENGTH 65
#endif
#ifndef SNAPERZ_PERIOD
#define SNAPERZ_PERIOD 12
#endif
#ifndef SNAPERZ_HARD_PUSH_LIMIT
#define SNAPERZ_HARD_PUSH_LIMIT 12
#endif

static constexpr uint32_t kLength = SNAPERZ_LENGTH;
static constexpr uint32_t kPeriod = SNAPERZ_PERIOD;
static constexpr uint32_t kHardPushLimit = SNAPERZ_HARD_PUSH_LIMIT;

// Constants. These are mirrored as macros because the choice of extender
// implementation is made by the preprocessor, which cannot see constexpr.
#define SNAPERZ_VIRTUAL_PUSH_LIMIT (SNAPERZ_PERIOD / 4 - 2)
#if SNAPERZ_VIRTUAL_PUSH_LIMIT < SNAPERZ_HARD_PUSH_LIMIT
#define SNAPERZ_PUSH_LIMIT SNAPERZ_VIRTUAL_PUSH_LIMIT
#else
#define SNAPERZ_PUSH_LIMIT SNAPERZ_HARD_PUSH_LIMIT
#endif

static constexpr uint32_t kVirtualPushLimit = SNAPERZ_VIRTUAL_PUSH_LIMIT;
static constexpr uint32_t kPushLimit = SNAPERZ_PUSH_LIMIT;
static constexpr uint32_t kLastPushLimit = std::min(kPushLimit + 1, kHardPushLimit);

typedef smallest_fit<kLength + 1>::type len_t;

// Definitions for checking loops. Use 1 for on, 0 for off.
#ifndef CHECK_LOOP
#define CHECK_LOOP 1
#endif
// Can be up to 2 times faster at finding loops, but slows down simulation slightly.
#ifndef FAST_LOOP_DETECTION
#define FAST_LOOP_DETECTION 1
#endif

// Definitions for logging status updates
#ifndef LOG_STATUS_UPDATES
#define LOG_STATUS_UPDATES 1
#endif
// Interval in number of pulses
#define LOGGING_INTERVAL UINT64_C(100000000)
