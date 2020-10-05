cmake_minimum_required(VERSION 2.8.8)

set(RES_FILES
   ${RES_FILES}
   ${CMAKE_CURRENT_SOURCE_DIR}/history.comp
   ${CMAKE_CURRENT_SOURCE_DIR}/volume.comp
   ${CMAKE_CURRENT_SOURCE_DIR}/blur_z.comp
   ${CMAKE_CURRENT_SOURCE_DIR}/blur_y.comp
   ${CMAKE_CURRENT_SOURCE_DIR}/blur_x.comp
   PARENT_SCOPE
)
