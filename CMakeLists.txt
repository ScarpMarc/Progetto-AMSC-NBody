# CMakeList.txt: file di progetto CMake di primo livello. Eseguire la configurazione globale
# e includere qui i sottoprogetti.
#
cmake_minimum_required (VERSION 3.20)

project ("Progetto AMSC NBody")

find_package(OpenGL REQUIRED)

# Includere i sottoprogetti.
add_subdirectory ("Progetto AMSC NBody")

############
# GRAPHICS #
############

#GLM
add_subdirectory("glm")
target_link_libraries(main glm::glm)

#GLFW
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
add_subdirectory("glfw")
target_link_libraries(main glfw)

#GLEW
add_subdirectory("glew/build/cmake")
target_link_libraries(main glew)

target_link_libraries(main OpenGL::GL)
