cmake_minimum_required(VERSION 3.20)

project("Progetto AMSC Nbody")

set(sources main.cpp Profiling.cpp)

if(OPENGL_FOUND)
set(sources ${sources} Controls.cpp OpenGLFunctions.cpp shader.cpp)
endif()

add_executable(main ${sources})

target_include_directories(main PUBLIC .)
set_property(TARGET main PROPERTY CXX_STANDARD 17)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")