# CMakeList.txt: file di progetto CMake di primo livello. Eseguire la configurazione globale
# e includere qui i sottoprogetti.
#
cmake_minimum_required (VERSION 3.20)

project ("Progetto AMSC NBody")

find_package(OpenGL OPTIONAL_COMPONENTS OpenGL)
if(OPENGL_FOUND)
add_compile_definitions(USE_GRAPHICS)
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