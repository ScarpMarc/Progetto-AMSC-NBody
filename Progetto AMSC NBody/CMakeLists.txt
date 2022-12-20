cmake_minimum_required(VERSION 3.20)

project("Progetto AMSC Nbody")

add_executable(main main.cpp)

target_include_directories(main PUBLIC .)
set_property(TARGET main PROPERTY CXX_STANDARD 17)
