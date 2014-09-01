SET(PROJECT_NAME globOpt)
PROJECT(${PROJECT_NAME})
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

SET( QCQPCPP_DIR "/home/bontius/workspace/QCQPcpp/" )
SET( MOSEK_DIR "/home/bontius/workspace/3rdparty/mosek/7/tools/platform/linux64x86/" )
SET( BONMIN_DIR "/home/bontius/workspace/3rdparty/CoinBonmin-stable/build/")
SET( GUROBI_DIR "/home/bontius/workspace/3rdparty/gurobi562/linux64/" )
#SET( FH_DIR /home/bontius/workspace/3rdparty/segmentFH )
SET( OpenCV_DIR "/home/bontius/workspace/3rdparty/opencv-trunk/install/share/OpenCV/" )
SET( PCL_DIR "/home/bontius/workspace/3rdparty/pcl-trunk/install/share/pcl-1.7/")
SET( PCLTOOLS_DIR "/home/bontius/workspace/geometry-tools/PCLTools/" )

# BUILD
SET( dir ${CMAKE_CURRENT_SOURCE_DIR}/build/${CMAKE_BUILD_TYPE} )
SET( EXECUTABLE_OUTPUT_PATH ${dir}/bin )#CACHE PATH "Build directory" FORCE)
SET( LIBRARY_OUTPUT_PATH ${dir}/lib )#CACHE PATH "Build directory" FORCE)

SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -Wall -fno-strict-aliasing -Wno-unused-local-typedefs" ) #-Wreturn-type
SET( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -g" )
SET( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O0" )
MESSAGE(STATUS "CXX_FLAGS: " ${CMAKE_CXX_FLAGS} "\nCXX_FLAGS" ${CMAKE_CXX_FLAGS_DEBUG} )

#________________DEPS________________#
#____________________________________#

#OpenGL
FIND_PACKAGE( OpenGL REQUIRED )
INCLUDE_DIRECTORIES( ${OPENGL_INCLUDE_DIRS} )

#OpenCV
FIND_PACKAGE( OpenCV COMPONENTS imgproc highgui core REQUIRED)
INCLUDE_DIRECTORIES( ${OpenCV_INCLUDE_DIRS} )

#PCL
FIND_PACKAGE( PCL REQUIRED ) # COMPONENTS COMMON IO SEARCH FEATURES FILTERS SEGMENTATION VISUALIZATION SAMPLE_CONSENSUS
INCLUDE_DIRECTORIES( ${PCL_INCLUDE_DIRS} )
ADD_DEFINITIONS( ${PCL_DEFINITIONS} )

# GUROBI - gurobi license in symlink /opt/gurobi/gurobi.lic
add_library( libgurobi SHARED IMPORTED )
set_property( TARGET libgurobi PROPERTY IMPORTED_LOCATION_RELEASE ${GUROBI_DIR}/lib/libgurobi56.so )
set_property( TARGET libgurobi PROPERTY IMPORTED_LOCATION_DEBUG ${GUROBI_DIR}/lib/libgurobi56.so )
add_library( libgurobi_c++ STATIC IMPORTED )
set_property( TARGET libgurobi_c++ PROPERTY IMPORTED_LOCATION_RELEASE ${GUROBI_DIR}/lib/libgurobi_c++.a )
set_property( TARGET libgurobi_c++ PROPERTY IMPORTED_LOCATION_DEBUG ${GUROBI_DIR}/lib/libgurobi_c++.a )
INCLUDE_DIRECTORIES( ${GUROBI_DIR}/include )
SET( GUROBI_LIBRARIES libgurobi_c++ libgurobi )

# Mosek
add_library( libmosek64 SHARED IMPORTED )
set_property( TARGET libmosek64 PROPERTY IMPORTED_LOCATION_RELEASE ${MOSEK_DIR}/bin/libmosek64.so )
set_property( TARGET libmosek64 PROPERTY IMPORTED_LOCATION_DEBUG ${MOSEK_DIR}/bin/libmosek64.so )
INCLUDE_DIRECTORIES( ${MOSEK_DIR}/h )
SET( MOSEK_LIBRARIES libmosek64 )

# Bonmin
SET( BONMIN_LIBRARIES
    ${BONMIN_DIR}/lib/libipopt.so
    ${BONMIN_DIR}lib/libbonmin.so
    ${BONMIN_DIR}lib/libCbc.so
    #${BONMIN_DIR}lib/libCbcSolver.so
    ${BONMIN_DIR}lib/libCgl.so
    ${BONMIN_DIR}lib/libClp.so
    ${BONMIN_DIR}lib/libcoinhsl.so
    ${BONMIN_DIR}lib/libCoinUtils.so
    ${BONMIN_DIR}lib/libOsi.so
    #${BONMIN_DIR}lib/libOsiCbc.so
    ${BONMIN_DIR}lib/libOsiClp.so
    #${BONMIN_DIR}lib/libOsiCommonTests.so
    lapack
    dl
)
INCLUDE_DIRECTORIES( ${BONMIN_DIR}/include )

# QCQPCPP
INCLUDE_DIRECTORIES( ${QCQPCPP_DIR}/include )

# libPCLTools
add_library( libPCLTools SHARED IMPORTED )
set_property( TARGET libPCLTools PROPERTY IMPORTED_LOCATION_RELEASE ${PCLTOOLS_DIR}/build/${CMAKE_BUILD_TYPE}/libPCLTools.so )
set_property( TARGET libPCLTools PROPERTY IMPORTED_LOCATION_DEBUG   ${PCLTOOLS_DIR}/build/${CMAKE_BUILD_TYPE}/libPCLToolsd.so )
INCLUDE_DIRECTORIES( ${PCLTOOLS_DIR}/include )
SET( PCLTOOLS_LIBRARIES libPCLTools )

#INCLUDE_DIRECTORIES( ${FH_DIR}/include )

#________________GLOBOPT________________#
#_______________________________________#

ADD_DEFINITIONS( -DGF2_USE_PCL=1
#                 -D__x86_64__=1
#                 -Dlinux=1
#                 -DHAVE_OPENCV=1
#                 #-DHAVE_OPENNI=1
#                 -DQT_NO_KEYWORDS
#                 -DVTK_EXCLUDE_STRSTREAM_HEADERS
                 )

SET( GLOBOPT_TARGET_NAME gurobi_opt )
INCLUDE_DIRECTORIES( include
    /usr/include/eigen3/
)

SET( GLOBOPT_H_LIST
    include/globfit2/my_types.h
    include/globfit2/util/diskUtil.hpp
    include/globfit2/util/util.hpp
    include/globfit2/util/pcl_util.hpp
    include/globfit2/primitives/linePrimitive2.h
    include/globfit2/primitives/pointPrimitive.h
    include/globfit2/visualization/visualizer.h
    include/globfit2/ground_truth/gtCreator.h
    include/globfit2/optimization/energyFunctors.h
    include/globfit2/optimization/candidateGenerator.h
    include/globfit2/optimization/merging.hpp
    include/optimization/qp/gurobiOpt.h
    include/optimization/qp/MyGRBCallback.h
    include/optimization/qp/solver.h
    include/globfit2/optimization/problemSetup.hpp
    include/params.h
    include/globfit2/optimization/patchDistanceFunctors.h
)

SET( GLOBOPT_SRC_LIST
    src/ground_truth/ground_truth.cpp
    src/ground_truth/gtCreator.cpp
    src/optimization/qp/qp_main.cpp
    src/optimization/qp/candidateGenerator.cpp # cached compilation
    src/optimization/qp/visualizer.cpp         # cached compilation
)

ADD_EXECUTABLE( ${GLOBOPT_TARGET_NAME} ${GLOBOPT_SRC_LIST} ${GLOBOPT_H_LIST} )
target_link_libraries( ${GLOBOPT_TARGET_NAME}
    ${GUROBI_LIBRARIES}
    ${MOSEK_LIBRARIES}
    ${BONMIN_LIBRARIES}
    ${PCLTOOLS_LIBRARIES}
    ${OpenCV_LIBRARIES}
    ${PCL_LIBRARIES}
    boost_filesystem
    boost_system
)

