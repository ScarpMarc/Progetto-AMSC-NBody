cmake_minimum_required(VERSION 3.20)

project("Progetto AMSC Nbody")

add_executable(main main.cpp Profiling.cpp Controls.cpp OpenGLFunctions.cpp shader.cpp)

target_include_directories(main PUBLIC .)
set_property(TARGET main PROPERTY CXX_STANDARD 17)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(main PUBLIC OpenMP::OpenMP_CXX)
endif()
