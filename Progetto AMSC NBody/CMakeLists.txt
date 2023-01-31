cmake_minimum_required(VERSION 3.20)

project("Progetto AMSC Nbody")

set(sources main.cpp Profiling.cpp)

if(OPENGL_FOUND)
set(sources ${sources} Controls.cpp OpenGLFunctions.cpp shader.cpp)
endif()

add_executable(NBody ${sources})

target_include_directories(NBody PUBLIC .)
set_property(TARGET NBody PROPERTY CXX_STANDARD 20)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")