cmake_minimum_required(VERSION 2.8.8)

set(RES_FILES
   ${RES_FILES}
   ${CMAKE_CURRENT_SOURCE_DIR}/voxel.vert
   ${CMAKE_CURRENT_SOURCE_DIR}/voxel.geom
   ${CMAKE_CURRENT_SOURCE_DIR}/voxel.frag
   PARENT_SCOPE
)
