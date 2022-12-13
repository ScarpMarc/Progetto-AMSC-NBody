#pragma once

///
///
/// COMPILER FLAGS
///
///

/// ADVANCED_PROFILING: more profiling info, slower execution (for testing purposes).
#define ADVANCED_PROFILING

#ifdef ADVANCED_PROFILING
#include <chrono>
#include <array>
#include "Constants.h"

long long int forceComp_mean_durations_per_tick = 0, posComp_mean_durations_per_tick = 0, matrixComp_mean_duration = 0;

#endif

#include <ctime>
time_t programme_start;