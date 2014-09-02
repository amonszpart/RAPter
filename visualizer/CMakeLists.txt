SET(PROJECT_NAME globOptVis)
PROJECT(${PROJECT_NAME})
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

SET( PCL_DIR "/home/bontius/workspace/3rdparty/pcl-trunk2/install/share/pcl-1.7/"
          CACHE FILEPATH "Folder containing \"PCLConfig.cmake\", usually \"{PCL_INSTALL_DIR}/share/pcl-1.7/\".")
SET( OpenCV_DIR "/home/bontius/workspace/3rdparty/opencv-trunk/install/share/OpenCV/"
        CACHE FILEPATH "Folder containing the OpenCVConfig.cmake file, usually \"{OPENCV_INSTALL_DIR}/share/OpenCV/\"" FORCE)
#SET( PATH_PCLTOOLS_ROOT "/home/bontius/workspace/geometry-tools/PCLTools/"
#        CACHE FILEPATH "Folder containing \"include\" directory from the SG PCLTools library. Can be downloaded from \"https://github.com/smartgeometry-ucl/geometry-tools/tree/master/PCLTools\".")

# BUILD
SET( dir ${CMAKE_CURRENT_SOURCE_DIR}/build/${CMAKE_BUILD_TYPE} )
SET( EXECUTABLE_OUTPUT_PATH ${dir}/bin CACHE PATH "Build directory" FORCE)
SET( LIBRARY_OUTPUT_PATH ${dir}/lib CACHE PATH "Build directory" FORCE)

#OpenGL
FIND_PACKAGE( OpenGL REQUIRED )
INCLUDE_DIRECTORIES( ${OPENGL_INCLUDE_DIRS} )

#PCL
FIND_PACKAGE( PCL REQUIRED ) # COMPONENTS COMMON IO SEARCH FEATURES FILTERS SEGMENTATION VISUALIZATION SAMPLE_CONSENSUS
INCLUDE_DIRECTORIES( ${PCL_INCLUDE_DIRS} )
ADD_DEFINITIONS( ${PCL_DEFINITIONS} )

ADD_DEFINITIONS( -DGF2_USE_PCL=1
                 -DVTK_EXCLUDE_STRSTREAM_HEADERS
                 )

#OpenCV
FIND_PACKAGE( OpenCV COMPONENTS imgproc highgui core REQUIRED)
INCLUDE_DIRECTORIES( ${OpenCV_INCLUDE_DIRS} )

include_directories(
    include
    ../globOpt/include
)

SET( GLOBOPTVIS_H_LIST
    include/globfit2/visualization/visualization.h
)

SET( GLOBOPTVIS_HPP_LIST
    include/globfit2/visualization/visualizer.hpp
)

SET( GLOBOPTVIS_SRC_LIST
    src/visualization.cpp
)

SET( LIB_NAME lib${PROJECT_NAME})
ADD_LIBRARY( ${LIB_NAME} STATIC ${GLOBOPTVIS_SRC_LIST} ${GLOBOPTVIS_HPP_LIST} ${GLOBOPTVIS_H_LIST})
SET_TARGET_PROPERTIES( ${LIB_NAME} PROPERTIES DEBUG_POSTFIX "d")
TARGET_LINK_LIBRARIES( ${LIB_NAME}
    ${PCL_LIBRARIES}
    ${OpenCV_LIBRARIES}
    boost_system
)

ADD_EXECUTABLE( ${PROJECT_NAME}
    ${GLOBOPTVIS_SRC_LIST}
    ${GLOBOPTVIS_HPP_LIST}
    src/main.cpp
    ${GLOBOPTVIS_H_LIST}
)
TARGET_LINK_LIBRARIES( ${PROJECT_NAME}
    ${PCL_LIBRARIES}
    ${OpenCV_LIBRARIES}
    boost_system
)

