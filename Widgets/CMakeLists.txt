cmake_minimum_required(VERSION 2.8.8)

set(SOURCE
   ${SOURCE}
   ${CMAKE_CURRENT_SOURCE_DIR}/mainwindow.cpp
   ${CMAKE_CURRENT_SOURCE_DIR}/glcanvas.cpp
   PARENT_SCOPE
)

set(HEADERS
   ${HEADERS}
   ${CMAKE_CURRENT_SOURCE_DIR}/mainwindow.h
   ${CMAKE_CURRENT_SOURCE_DIR}/glcanvas.h
   PARENT_SCOPE
)

set(UI
    ${UI}
    ${CMAKE_CURRENT_SOURCE_DIR}/mainwindow.ui
    PARENT_SCOPE
)
