cmake_minimum_required(VERSION 2.8.8)

set(SOURCE
   ${SOURCE}
   ${CMAKE_CURRENT_SOURCE_DIR}/psfsolver.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/psfsolvergpu.cpp
   PARENT_SCOPE
)

set(HEADERS
   ${HEADERS}
   ${CMAKE_CURRENT_SOURCE_DIR}/psfsolver.h
   ${CMAKE_CURRENT_SOURCE_DIR}/psfsolvergpu.h
   PARENT_SCOPE
)
