#pragma once

#include <vector>

/// <summary>
/// Simple struct that defines a colour
/// </summary>
struct Colour
{

	unsigned char r, g, b, a;
};

/// <summary>
/// Defines visualisation info for particles and derived types. 
/// </summary>
struct VisualisationInfo
{
	/// <summary>
	/// Colour palette. If the amount of colour is greater than 1, the colours will cycle.
	/// </summary>
	std::vector<Colour> colours;

	/// <summary>
	/// Size for OpenGL rendering.
	/// </summary>
	float size;

	/// <summary>
	/// Amount of frames to wait for when cycling colours.
	/// </summary>
	unsigned int colour_cycle_wait_time;

	/// <summary>
	/// Whether the particle should be rendered at all.
	/// </summary>
	bool invisible;
};

/// <summary>
/// Defines a list of VisualisationInfos that should go together. 
/// Used upon colour palette change to assign random effects to particles.
/// For example, one may want to view the particles as stars, which are from yellow to blue and glowing;
///		the initialisation function will take care to assign to each particle a different, random VisualisationInfo.
/// </summary>
struct VisualisationInfo_palette
{
	std::vector<VisualisationInfo> choices;
};

// DEFINITIONS 

/*
 From http://www.vendian.org/mncharity/dir3/starcolor/
	O5(V)  	  157 180 255   #9db4ff
	B1(V)  	  162 185 255   #a2b9ff
	B3(V)  	  167 188 255   #a7bcff
	B5(V)  	  170 191 255   #aabfff
	B8(V)  	  175 195 255   #afc3ff
	A1(V)  	  186 204 255   #baccff
	A3(V)  	  192 209 255   #c0d1ff
	A5(V)  	  202 216 255   #cad8ff
	F0(V)  	  228 232 255   #e4e8ff
	F2(V)  	  237 238 255   #edeeff
	F5(V)  	  251 248 255   #fbf8ff
	F8(V)  	  255 249 249   #fff9f9
	G2(V)  	  255 245 236   #fff5ec
	G5(V)  	  255 244 232   #fff4e8
	G8(V)  	  255 241 223   #fff1df
	K0(V)  	  255 235 209   #ffebd1
	K4(V)  	  255 215 174   #ffd7ae
	K7(V)  	  255 198 144   #ffc690
	M2(V)  	  255 190 127   #ffbe7f
	M4(V)  	  255 187 123   #ffbb7b
	M6(V)  	  255 187 123   #ffbb7b
*/