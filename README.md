# Progetto-AMSC-NBody

# How to run
## Mandatory options
The user must select at one of the two
- `-json` to load settings from `settings.json` file and generate particles randomly
- `-load filename` to load particles from a file

## Required files
`vertexShader.vertexshader`, `fragmentShader.fragmentshader` and `particle.dds` must be in the same directory as the executable
Optionally, the executable must have `settings.json` in the same directory if the user wants a JSON generation.

## Controls
Left and right arrows to strafe left and right; up and down arrows to zoom; ESC to close the graphics window. This will halt the entire programme.

# How to build 
Skip sections 1 and 2 if you _do not_ want graphics.

**WARNING**: if you have installed the [mk-modules](https://github.com/elauksap/mk) and want to use the graphics, you **MUST** unload `gcc-glibc` as this interferes with the CMake routines that find OpenGL. In general, you should always refer to the latest version of CMake.
1. Go to `glew/auto` and `make`; this will configure the appropriate dependencies. You may need to access `Makefile` and `modify PYTHON ?= python` into `PYTHON ?= python3`
2. You may need to install the following libraries: `libxrandr-dev libxinerama-dev libxcursor-dev libxi-dev`
3. Create a `build` directory in the root folder and `cmake ..` inside of it. The output will state whether OpenGL has been found (if not, read the warning and refer to points 1 and 2).
4. `make`. The result should be in the `Progetto AMSC NBody` folder in `build`
