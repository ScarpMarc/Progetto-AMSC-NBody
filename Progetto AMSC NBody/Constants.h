#pragma once

//
// MISCELLANEOUS
//

/// <summary>
/// Particle dimension (2D, 3D)
/// </summary>
const unsigned int DIM = 3;

const unsigned int total_particles = 20;

//
// TIMING
//

/// <summary>
/// Amount of ticks computed during the simulation
/// </summary>
const unsigned long int max_ticks = 100000;

/// <summary>
/// Amount of tick per simulation second.
/// Note: this has no connection to real-time seconds, but is only simulation-related.
/// </summary>
const double ticks_per_second = 100;

//
// GRAPHICS
//

/// <summary>
/// Amount of real-time milliseconds that the drawing routine pauses for after each drawing cycle.
/// </summary>
const unsigned int screen_refresh_millis = 200;

/// <summary>
/// Screen resolution in the horizontal direction
/// </summary>
const unsigned int screenResX = 1024;
/// <summary>
/// Screen resolution in the vertical direction
/// </summary>
const unsigned int screenResY = 1024;