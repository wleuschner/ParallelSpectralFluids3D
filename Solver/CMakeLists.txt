cmake_minimum_required(VERSION 2.8.8)

add_subdirectory("SolverImpl")

set(SOURCE
   ${SOURCE}
   ${CMAKE_CURRENT_SOURCE_DIR}/abstractsolver.cpp
   PARENT_SCOPE
)

set(HEADERS
   ${HEADERS}
   ${CMAKE_CURRENT_SOURCE_DIR}/abstractsolver.h
   PARENT_SCOPE
)
