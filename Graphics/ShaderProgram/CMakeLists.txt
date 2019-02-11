cmake_minimum_required(VERSION 2.8.8)

set(SOURCE
   ${SOURCE}
   ${CMAKE_CURRENT_SOURCE_DIR}/ShaderProgram.cpp
   PARENT_SCOPE
)

set(HEADERS
   ${HEADERS}
   ${CMAKE_CURRENT_SOURCE_DIR}/ShaderProgram.h
   PARENT_SCOPE
)
