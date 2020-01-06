cmake_minimum_required(VERSION 2.8.8)

add_subdirectory("Phong")
add_subdirectory("Line")
add_subdirectory("Points")
add_subdirectory("Volume")

set(RES_FILES
   ${RES_FILES}
   PARENT_SCOPE
)
