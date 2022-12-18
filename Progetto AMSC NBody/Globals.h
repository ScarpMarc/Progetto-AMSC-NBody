#pragma once

///
///
/// COMPILER FLAGS
///
///

/// ADVANCED_PROFILING: more profiling info, slower execution (for testing purposes).
#define ADVANCED_PROFILING

#ifdef ADVANCED_PROFILING
//#include <chrono>
//#include <array>
#include "Constants.h"

// All in microseconds
extern long long int forceComp_mean_durations_per_tick, posComp_mean_durations_per_tick, matrixComp_mean_duration;
#endif

extern long long int total_sim_duration;


#include <ctime>
extern time_t programme_start;

