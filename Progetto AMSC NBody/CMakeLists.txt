cmake_minimum_required(VERSION 3.20)

project("Progetto AMSC Nbody")

add_executable(main main.cpp Profiling.cpp Controls.cpp OpenGLFunctions.cpp shader.cpp)

target_include_directories(main PUBLIC .)
set_property(TARGET main PROPERTY CXX_STANDARD 17)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")