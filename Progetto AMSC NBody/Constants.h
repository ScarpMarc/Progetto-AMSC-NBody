#pragma once

//
// MISCELLANEOUS
//

/// <summary>
/// Particle dimension (2D, 3D)
/// </summary>
const unsigned int DIM = 3;

extern unsigned int total_particles;

//
// TIMING
//

/// <summary>
/// Amount of ticks computed during the simulation
/// </summary>
extern unsigned long int max_ticks;

/// <summary>
/// Amount of tick per simulation second.
/// Note: this has no connection to real-time seconds, but is only simulation-related.
/// </summary>
extern double ticks_per_second;

/// <summary>
/// Amount of steps without saving particles status.
/// Note: this has no connection to real-time seconds, but is only simulation-related.
/// </summary>
extern unsigned int save_status_interval;

//
// GRAPHICS
//

/// <summary>
/// Amount of real-time milliseconds that the drawing routine pauses for after each drawing cycle.
/// </summary>
extern unsigned int screen_refresh_millis;

/// <summary>
/// Screen resolution in the horizontal direction
/// </summary>
extern unsigned int screenResX;
/// <summary>
/// Screen resolution in the vertical direction
/// </summary>
extern unsigned int screenResY;