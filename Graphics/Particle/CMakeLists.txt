cmake_minimum_required(VERSION 2.8.8)

set(HEADERS
   ${HEADERS}
   ${CMAKE_CURRENT_SOURCE_DIR}/Particle.cpp
   PARENT_SCOPE
)

set(SOURCE
   ${SOURCE}
   ${CMAKE_CURRENT_SOURCE_DIR}/Particle.h
   PARENT_SCOPE
)
