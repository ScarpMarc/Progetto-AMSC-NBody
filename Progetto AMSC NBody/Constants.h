#pragma once

//
// MISCELLANEOUS
//

/// <summary>
/// Particle dimension (2D, 3D)
/// </summary>
const unsigned int DIM = 3;

/// <summary>
/// Particle dimension (2D, 3D)
/// </summary>
unsigned int total_particles = 3;

//
// TIMING
//

/// <summary>
/// Amount of ticks computed during the simulation
/// </summary>
unsigned long int max_ticks = 10000;

/// <summary>
/// Amount of tick per simulation second.
/// Note: this has no connection to real-time seconds, but is only simulation-related.
/// </summary>
double ticks_per_second = 10;

/// <summary>
/// Amount of steps without saving particles status.
/// Note: this has no connection to real-time seconds, but is only simulation-related.
/// </summary>
unsigned int save_status_interval = 10;

//
// GRAPHICS
//

/// <summary>
/// Amount of real-time milliseconds that the drawing routine pauses for after each drawing cycle.
/// </summary>
unsigned int screen_refresh_millis = 100;

/// <summary>
/// Screen resolution in the horizontal direction
/// </summary>
unsigned int screenResX = 1024;
/// <summary>
/// Screen resolution in the vertical direction
/// </summary>
unsigned int screenResY = 1024;
