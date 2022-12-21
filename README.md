# Progetto-AMSC-NBody

# How to run
## Mandatory options
The user must select at one of the two
- `-json` to load settings from `settings.json` file and generate particles randomly
- `-load filename` to load particles from a file

## Graphics
Goto `glew/auto` and make; this will configure the appropriate dependencies.

## Required files
`vertexshader.vertexshader`, `fragmentshader.fragmentshader` and `particle.dds` must be in the same directory as the executable
Optionally, the executable must have `settings.json` in the same directory if the user wants a JSON generation.

# Controls
Left and right arrows to strafe left and right; up and down arrows to zoom; ESC to close the graphics window. This will halt the entire programme.
