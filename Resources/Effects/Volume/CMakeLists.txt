cmake_minimum_required(VERSION 2.8.8)

set(RES_FILES
   ${RES_FILES}
   ${CMAKE_CURRENT_SOURCE_DIR}/volume.vert
   ${CMAKE_CURRENT_SOURCE_DIR}/volume.frag
   ${CMAKE_CURRENT_SOURCE_DIR}/blit.frag
   PARENT_SCOPE
)
