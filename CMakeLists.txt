# CMakeList.txt: file di progetto CMake di primo livello. Eseguire la configurazione globale
# e includere qui i sottoprogetti.
#
cmake_minimum_required (VERSION 3.20)
include(${CMAKE_CURRENT_SOURCE_DIR}/CMakeVarious/OptimizeForArchitecture.cmake)

project ("Progetto AMSC NBody")

OptimizeForArchitecture()

SET(NO_GRAPHICS "OFF" CACHE BOOL "Force graphics OFF")

if(NO_GRAPHICS)
    message("Compiling WITHOUT graphics")
endif()

option(SCOREP_ENABLE "Create PIConGPU with Score-P support" OFF)

if(SCOREP_ENABLE)
    message(STATUS "Building with Score-P support")
    set(SCOREP_ROOT "$ENV{SCOREP_ROOT}")
    if(NOT SCOREP_ROOT)
        message(FATAL_ERROR "Environment variable SCOREP_ROOT not set!")
    endif(NOT SCOREP_ROOT)

    # compile flags
    execute_process(COMMAND $ENV{SCOREP_ROOT}/bin/scorep-config --nocompiler --cflags
                    OUTPUT_VARIABLE SCOREP_COMPILEFLAGS
                    RESULT_VARIABLE SCOREP_CONFIG_RETURN)
    if(NOT SCOREP_CONFIG_RETURN EQUAL 0)
        message(FATAL_ERROR "Can NOT execute 'scorep-config' at $ENV{SCOREP_ROOT}/bin/scorep-config - check file permissions")
    endif()

    # link flags
    execute_process(COMMAND $ENV{SCOREP_ROOT}/bin/scorep-config --mpp=none --thread=omp --ldflags
                    OUTPUT_VARIABLE SCOREP_LINKFLAGS)
    # libraries
    execute_process(COMMAND $ENV{SCOREP_ROOT}/bin/scorep-config --mpp=none --thread=omp --libs
                    OUTPUT_VARIABLE SCOREP_LIBFLAGS)
    string(STRIP "${SCOREP_LIBFLAGS}" SCOREP_LIBFLAGS)

    # subsystem iniialization file
    execute_process(COMMAND $ENV{SCOREP_ROOT}/bin/scorep-config --mpp=none --thread=omp --adapter-init
                    OUTPUT_VARIABLE SCOREP_INIT_FILE)
    file(WRITE ${CMAKE_BINARY_DIR}/scorep_init.c "${SCOREP_INIT_FILE}")

    if(SCOREP_ENABLE)
        set(SCOREP_SRCFILES "${CMAKE_BINARY_DIR}/scorep_init.c")
    endif(SCOREP_ENABLE)

    # modify our flags
    set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${SCOREP_LINKFLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SCOREP_COMPILEFLAGS}")
endif(SCOREP_ENABLE)

if(NOT NO_GRAPHICS)
    find_package(OpenGL OPTIONAL_COMPONENTS OpenGL)
    if(OPENGL_FOUND)
        add_compile_definitions(USE_GRAPHICS)
    endif()
endif()
#target_link_libraries(NBody PUBLIC OpenGL::OpenGL)

# Includere i sottoprogetti.
add_subdirectory ("Progetto AMSC NBody")

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(NBody PUBLIC OpenMP::OpenMP_CXX)
endif()

############
# GRAPHICS #
############

if(OPENGL_FOUND)
if(TARGET OpenGL::OpenGL)
    target_link_libraries(NBody PUBLIC OpenGL::OpenGL)
else()
    target_link_libraries(NBody PUBLIC OpenGL::GL)
endif()

#GLM
add_subdirectory("glm")
target_link_libraries(NBody PUBLIC glm::glm)

#GLFW
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
add_subdirectory("glfw")
target_link_libraries(NBody PUBLIC glfw)

#GLEW
add_subdirectory("glew/build/cmake")
target_link_libraries(NBody PUBLIC glew)
endif()

add_custom_target("deploy")
add_custom_command(TARGET "deploy" POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy_directory
${CMAKE_CURRENT_SOURCE_DIR}/Resources
${CMAKE_CURRENT_BINARY_DIR}/'Progetto AMSC NBody')