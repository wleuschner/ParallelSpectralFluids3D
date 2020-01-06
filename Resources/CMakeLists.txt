cmake_minimum_required(VERSION 2.8.8)

add_subdirectory("Effects")
add_subdirectory("Voxel")
add_subdirectory("Simulation")
add_subdirectory("Volume")

set(RES_FILES
   ${RES_FILES}
   PARENT_SCOPE
)
