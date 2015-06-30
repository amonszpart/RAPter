SET(PROJECT_NAME RAPterVis)
PROJECT(${PROJECT_NAME})
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#_____________________________________#
#_______________PATHS_________________#
#_____________________________________#

SET( WORKSPACE_DIR "$ENV{HOME}/workspace" )
SET( THIRD_PARTY_DIR "${WORKSPACE_DIR}/3rdparty" )

# PCL
SET( PCL_DIR "${THIRD_PARTY_DIR}/pcl-trunk2/install/share/pcl-1.8/"
          CACHE FILEPATH "Folder containing \"PCLConfig.cmake\", usually \"{PCL_INSTALL_DIR}/share/pcl-1.8/\".")

#_____________________________________#
#________________BUILD________________#
#_____________________________________#

SET( dir ${CMAKE_CURRENT_SOURCE_DIR}/build/${CMAKE_BUILD_TYPE} )
SET( EXECUTABLE_OUTPUT_PATH ${dir}/bin CACHE PATH "Build directory" FORCE)
SET( LIBRARY_OUTPUT_PATH ${dir}/lib CACHE PATH "Build directory" FORCE)

SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-function -std=c++0x" )

#____________________________________#
#________________DEPS________________#
#____________________________________#

# OpenGL
FIND_PACKAGE( OpenGL REQUIRED )
INCLUDE_DIRECTORIES( ${OPENGL_INCLUDE_DIRS} )

#PCL
FIND_PACKAGE( PCL REQUIRED ) # COMPONENTS COMMON IO SEARCH FEATURES FILTERS SEGMENTATION VISUALIZATION SAMPLE_CONSENSUS
INCLUDE_DIRECTORIES( ${PCL_INCLUDE_DIRS} )
ADD_DEFINITIONS( ${PCL_DEFINITIONS} )

ADD_DEFINITIONS( -DRAPTER_USE_PCL=1
                 -DVTK_EXCLUDE_STRSTREAM_HEADERS
                 )

include_directories(
    include
    /usr/include/eigen3
    ../RAPter/include
)

SET( RAPTER_H_LIST
    include/rapter/visualization/visualization.h
)

SET( RAPTER_HPP_LIST
    include/rapter/visualization/visualizer.hpp
    include/rapter/visualization/impl/visualization.hpp
)

SET( RAPTER_SRC_LIST
    src/visualization.cpp
    ../RAPter/src/templateInstantiation/taggable.cpp
    ../RAPter/src/templateInstantiation/planePrimitive.cpp
    ../RAPter/src/templateInstantiation/linePrimitive.cpp
)

SET( LIB_NAME ${PROJECT_NAME} )
ADD_LIBRARY( ${LIB_NAME} STATIC ${RAPTER_SRC_LIST} ${RAPTER_HPP_LIST} ${RAPTER_H_LIST})
SET_TARGET_PROPERTIES( ${LIB_NAME} PROPERTIES DEBUG_POSTFIX "d")
TARGET_LINK_LIBRARIES( ${LIB_NAME}
                       ${PCL_LIBRARIES}
                       boost_system
)

SET( BIN_NAME run${PROJECT_NAME} )
ADD_EXECUTABLE( ${BIN_NAME}
                ${RAPTER_SRC_LIST}
                ${RAPTER_HPP_LIST}
                src/main.cpp
                ${RAPTER_H_LIST}
)

TARGET_LINK_LIBRARIES( ${BIN_NAME}
                       ${PCL_LIBRARIES}
                       boost_system
)

SET_TARGET_PROPERTIES( ${BIN_NAME}
                       PROPERTIES OUTPUT_NAME rapterVis )

