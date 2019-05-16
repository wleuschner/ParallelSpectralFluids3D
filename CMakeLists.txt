cmake_minimum_required(VERSION 2.8.8)

project(PSF3D)
set (CMAKE_CXX_STANDARD 11)
add_definitions(-DQT_NO_OPENGL_ES_2 -DGLM_FORCE_SSE2 -DGLM_FORCE_ALIGNED -DGLM_ENABLE_EXPERIMENTAL -DVIENNACL_HAVE_EIGEN -DMICROPROFILE_IMPL -DMICROPROFILE_WEBSERVER=1 -DMICROPROFILE_ENABLED=1 -DMICROPROFILE_GPU_TIMERS=0)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

#set(Qt5Core_DIR /opt/Qt/5.6.3/gcc_64/lib/cmake/Qt5Core)
#set(Qt5Gui_DIR /opt/Qt/5.6.3/gcc_64/lib/cmake/Qt5Gui)
#set(Qt5Widgets_DIR /opt/Qt/5.6.3/gcc_64/lib/cmake/Qt5Widgets)
#set(Qt5OpenGL_DIR /opt/Qt/5.6.3/gcc_64/lib/cmake/Qt5OpenGL)

find_package(Qt4 REQUIRED)
include(${QT_USE_FILE})
add_definitions( ${QT_DEFINITIONS} )
include_directories( ${CMAKE_BINARY_DIR} )
#find_package(Qt5Core)
#find_package(Qt5Gui)
#find_package(Qt5Widgets)
#find_package(Qt5OpenGL)

add_subdirectory("Solver")
add_subdirectory("DEC")
add_subdirectory("Geometry")
add_subdirectory("Graphics")
add_subdirectory("Widgets")
add_subdirectory("Resources")


set(CMAKE_INCLUDE_CURRENT_DIR ON)
#set(CMAKE_AUTOMOC ON)
#set(CMAKE_AUTOUIC ON)
#set(CMAKE_AUTORCC ON)

find_package(GLEW REQUIRED)
find_package(OpenGL REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(assimp)
find_package(OpenCL REQUIRED)

#include_directories(${Qt5Widgets_INCLUDES})
#add_definitions(${Qt5Widgets_DEFINITIONS})

set(UI
    ${UI}
)

set(HEADERS
   ${HEADERS}
   microprofile/microprofile.h
)


set(SOURCE
    ${SOURCE}
    main.cpp
    microprofile/microprofile.cpp
)

set(RES_FILES
    ${RES_FILES}
)

QT4_WRAP_UI(QT_HEADER_FORM Widgets/mainwindow.ui)
QT4_WRAP_CPP(QT_HEADER_MOC Widgets/mainwindow.h Widgets/glcanvas.h)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -mtune=native -march=native")

add_executable(${PROJECT_NAME} ${RES_FILES} ${UI} ${SOURCE} ${HEADERS} ${QT_HEADER_FORM} ${QT_HEADER_MOC})

#qt5_use_modules(${PROJECT_NAME} Widgets OpenGL Core Gui)

target_link_libraries(${PROJECT_NAME} ${QT_QTCORE_LIBRARY} ${QT_QTGUI_LIBRARY} ${QT_QTOPENGL_LIBRARY} ${QT_QTWIDGETS_LIBRARY} ${OPENGL_gl_LIBRARY} ${GLEW_LIBRARIES} ${EIGEN_LIBRARIES} ${OpenCL_LIBRARY} arpack assimp pthread)

set (source ${CMAKE_SOURCE_DIR}/Resources)
set (destination ${CMAKE_CURRENT_BINARY_DIR}/Res)
add_custom_command(
 TARGET ${PROJECT_NAME} POST_BUILD
 COMMAND ${CMAKE_COMMAND} -E create_symlink ${source} ${destination}
 DEPENDS ${destination}
 COMMENT "symbolic link resources folder from ${source} => ${destination}"
)
