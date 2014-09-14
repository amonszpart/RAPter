SET(PROJECT_NAME globOpt)
PROJECT(${PROJECT_NAME})
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

SET( WORKSPACE_DIR "$ENV{HOME}/workspace" )
SET( THIRD_PARTY_DIR "${WORKSPACE_DIR}/3rdparty" )

SET( OpenCV_DIR "${THIRD_PARTY_DIR}/opencv-trunk/install/share/OpenCV/"
        CACHE FILEPATH "Folder containing the OpenCVConfig.cmake file, usually \"{OPENCV_INSTALL_DIR}/share/OpenCV/\"")
SET( PCL_DIR "${THIRD_PARTY_DIR}/pcl-trunk2/install/share/pcl-1.7/"
          CACHE FILEPATH "Folder containing \"PCLConfig.cmake\", usually \"{PCL_INSTALL_DIR}/share/pcl-1.7/\".")

SET( PATH_MOSEK_DIR  "${THIRD_PARTY_DIR}/mosek/7/tools/platform/linux64x86" )
SET( PATH_BONMIN_DIR "${THIRD_PARTY_DIR}/CoinBonmin-stable/build/")

SET( PATH_GUROBI_DIR "${THIRD_PARTY_DIR}/gurobi562/linux64/"
        CACHE FILEPATH "Folder containing libgurobixx.so and libgurobi_c++.a, usually \"{GUROBI_INSTALL_DIR}/gurobi562/linux64/\".")

SET( PATH_LIBFBI_ROOT "${THIRD_PARTY_DIR}/libfbi/")

SET( PATH_QCQPCPP_ROOT "${WORKSPACE_DIR}/QCQPcpp/"
        CACHE FILEPATH "Root directory of QCQPcpp library containing the \"include\" directory.")
SET( PATH_PCLTOOLS_ROOT "${WORKSPACE_DIR}/geometry-tools/PCLTools/"
        CACHE FILEPATH "Folder containing \"include\" directory from the SG PCLTools library. Can be downloaded from \"https://github.com/smartgeometry-ucl/geometry-tools/tree/master/PCLTools\".")
SET( PATH_VIS_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../visualization/"
        CACHE FILEPATH "Folder containing separate visulizer library")

SET( WITH_MOSEK OFF CACHE BINARY "Allow Mosek to be used for optimization.")
IF(WITH_MOSEK)
    ADD_DEFINITIONS(-DGF2_WITH_MOSEK)
ENDIF()
SET( WITH_BONMIN ON CACHE BINARY "Allow Bonmin to be used for optimization.")
IF(WITH_BONMIN)
    ADD_DEFINITIONS(-DGF2_WITH_BONMIN)
ENDIF()
SET( WITH_GUROBI OFF CACHE BINARY "Allow Gruboi to be used for optimization.")
IF(WITH_GUROBI)
    ADD_DEFINITIONS(-DGF2_WITH_GUROBI)
ENDIF()

SET( WITH_SAMPLE_INPUT OFF CACHE BINARY "Compile --sample-input option")
IF(WITH_SAMPLE_INPUT)
    ADD_DEFINITIONS(-DGF2_WITH_SAMPLE_INPUT)
ENDIF(WITH_SAMPLE_INPUT)

SET( WITH_FULL_LINKAGE OFF CACHE BINARY "Compile other patchpatch distance functors than representativeSqr")
IF(WITH_FULL_LINKAGE)
    ADD_DEFINITIONS(-DGF2_WITH_FULL_LINKAGE)
ENDIF(WITH_FULL_LINKAGE)

#SET( FH_DIR ${THIRD_PARTY_DIR}/segmentFH )

# BUILD
SET( dir ${CMAKE_CURRENT_SOURCE_DIR}/build/${CMAKE_BUILD_TYPE} )
SET( EXECUTABLE_OUTPUT_PATH ${dir}/bin CACHE PATH "Build directory" FORCE)
SET( LIBRARY_OUTPUT_PATH ${dir}/lib CACHE PATH "Build directory" FORCE)

SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -Wall -Wno-unused-local-typedefs " ) #-Wreturn-type
#SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall" ) #-Wreturn-type
#SET( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG} " )
SET( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE} -O0" )

#________________DEPS________________#
#____________________________________#

#OpenGL
FIND_PACKAGE( OpenGL REQUIRED )
INCLUDE_DIRECTORIES( ${OPENGL_INCLUDE_DIRS} )

#OpenCV
FIND_PACKAGE( OpenCV COMPONENTS imgproc highgui core REQUIRED)
INCLUDE_DIRECTORIES( ${OpenCV_INCLUDE_DIRS} )

#PCL
if(1)
    FIND_PACKAGE( PCL REQUIRED ) # COMPONENTS COMMON IO SEARCH FEATURES FILTERS SEGMENTATION VISUALIZATION SAMPLE_CONSENSUS
    INCLUDE_DIRECTORIES( ${PCL_INCLUDE_DIRS} )
    ADD_DEFINITIONS( ${PCL_DEFINITIONS} )
else()
    SET( PCL_INCLUDE_DIRS
        /usr/include/eigen3/
        /usr/include/pcl-1.7/
        /usr/include/vtk-5.8/
    )
    SET( PCL_LIBRARIES
        pcl_visualization
        pcl_sample_consensus
        pcl_features
        pcl_segmentation
        pcl_common
        pcl_io
        pcl_search
        pcl_filters
        vtkCommon
        vtkHybrid
        vtkFiltering
        vtkRendering
        vtkGraphics
    )
endif()

INCLUDE_DIRECTORIES( ${PCL_INCLUDE_DIRS} )

IF(WITH_GUROBI)
    # GUROBI - gurobi license in symlink /opt/gurobi/gurobi.lic
    add_library( libgurobi SHARED IMPORTED )
    set_property( TARGET libgurobi PROPERTY IMPORTED_LOCATION_RELEASE ${PATH_GUROBI_DIR}/lib/libgurobi56.so )
    set_property( TARGET libgurobi PROPERTY IMPORTED_LOCATION_DEBUG ${PATH_GUROBI_DIR}/lib/libgurobi56.so )
    add_library( libgurobi_c++ STATIC IMPORTED )
    set_property( TARGET libgurobi_c++ PROPERTY IMPORTED_LOCATION_RELEASE ${PATH_GUROBI_DIR}/lib/libgurobi_c++.a )
    set_property( TARGET libgurobi_c++ PROPERTY IMPORTED_LOCATION_DEBUG ${PATH_GUROBI_DIR}/lib/libgurobi_c++.a )
    INCLUDE_DIRECTORIES( ${PATH_GUROBI_DIR}/include )
    SET( GUROBI_LIBRARIES libgurobi_c++ libgurobi )
ENDIF(WITH_GUROBI)

# Mosek
IF(WITH_MOSEK)
    add_library( libmosek64 SHARED IMPORTED )
    set_property( TARGET libmosek64 PROPERTY IMPORTED_LOCATION_RELEASE ${PATH_MOSEK_DIR}/bin/libmosek64.so )
    set_property( TARGET libmosek64 PROPERTY IMPORTED_LOCATION_DEBUG ${PATH_MOSEK_DIR}/bin/libmosek64.so )
    INCLUDE_DIRECTORIES( ${PATH_MOSEK_DIR}/h )
    SET( MOSEK_LIBRARIES "libmosek64" )
ENDIF(WITH_MOSEK)


# Bonmin
IF(WITH_BONMIN)
    IF(EXISTS ${PATH_BONMIN_DIR}lib/libcoinhsl.so)
      SET (BONMIN_SOLVER_LIB ${PATH_BONMIN_DIR}lib/libcoinhsl.so)
      MESSAGE(STATUS "HSL Solver available, use instead of MUMPS")
    ELSE()
        SET (BONMIN_SOLVER_LIB ${PATH_BONMIN_DIR}lib/libcoinmumps.so)
    ENDIF(EXISTS ${PATH_BONMIN_DIR}lib/libcoinhsl.so)

    SET( BONMIN_LIBRARIES
        ${PATH_BONMIN_DIR}/lib/libipopt.so
        ${PATH_BONMIN_DIR}/lib/libbonmin.so
        ${PATH_BONMIN_DIR}/lib/libCbc.so
        #${PATH_BONMIN_DIR}lib/libCbcSolver.so
        ${PATH_BONMIN_DIR}/lib/libCgl.so
        ${PATH_BONMIN_DIR}/lib/libClp.so
        ${BONMIN_SOLVER_LIB}
        #${BONMIN_DIR}lib/libcoinhsl.so
        ${PATH_BONMIN_DIR}/lib/libCoinUtils.so
        ${PATH_BONMIN_DIR}/lib/libOsi.so
        #${PATH_BONMIN_DIR}/lib/libOsiCbc.so
        ${PATH_BONMIN_DIR}/lib/libOsiClp.so
        #${PATH_BONMIN_DIR}/lib/libOsiCommonTests.so
        lapack
        dl
    )
    INCLUDE_DIRECTORIES( ${PATH_BONMIN_DIR}/include )
ENDIF(WITH_BONMIN)

# QCQPCPP
INCLUDE_DIRECTORIES( ${PATH_QCQPCPP_ROOT}/include )

# libPCLTools
INCLUDE_DIRECTORIES( ${PATH_PCLTOOLS_ROOT}/include )

# libfbi
INCLUDE_DIRECTORIES( ${PATH_LIBFBI_ROOT}/include )

# vis
IF(1)
    add_library( libglobOptVis STATIC IMPORTED )
    set_property( TARGET libglobOptVis PROPERTY IMPORTED_LOCATION_RELEASE "${PATH_VIS_ROOT}/build/Release/lib/libglobOptVis.a" )
    set_property( TARGET libglobOptVis PROPERTY IMPORTED_LOCATION_DEBUG ${PATH_VIS_ROOT}/build/Debug/lib/libglobOptVisd.a )
    INCLUDE_DIRECTORIES( ${PATH_VIS_ROOT}/include )
    SET( GLOBOPTVIS_LIBRARIES "libglobOptVis" )
ENDIF()


#INCLUDE_DIRECTORIES( ${FH_DIR}/include )

#________________GLOBOPT________________#
#_______________________________________#

ADD_DEFINITIONS( -DGF2_USE_PCL=1
#                 -D__x86_64__=1
#                 -Dlinux=1
#                 -DHAVE_OPENCV=1
#                 #-DHAVE_OPENNI=1
#                 -DQT_NO_KEYWORDS
                 -DVTK_EXCLUDE_STRSTREAM_HEADERS
                 )

SET( GLOBOPT_TARGET_NAME glob_opt )
INCLUDE_DIRECTORIES( include
    ${EIGEN_INCLUDE_DIRS}
)

SET( GLOBOPT_HPP_LIST
    include/globfit2/optimization/impl/problemSetup.hpp
    include/globfit2/optimization/impl/merging.hpp
    include/globfit2/util/diskUtil.hpp
    include/globfit2/util/util.hpp
    include/globfit2/util/pcl_util.hpp
    include/globfit2/util/containers.hpp
    include/globfit2/processing/util.hpp
)

SET( GLOBOPT_H_LIST
    #include/params.h
    include/globfit2/util/parse.h
    include/globfit2/parameters.h
    include/globfit2/my_types.h
    include/globfit2/globOpt_types.h
    include/globfit2/primitives/linePrimitive.h
    include/globfit2/primitives/primitive.h
    include/globfit2/primitives/linePrimitive2.h
    include/globfit2/primitives/pointPrimitive.h
    include/globfit2/primitives/planePrimitive.h
    include/globfit2/optimization/patchDistanceFunctors.h
    include/globfit2/optimization/energyFunctors.h
    include/globfit2/optimization/candidateGenerator.h
    include/globfit2/optimization/problemSetup.h
    include/globfit2/optimization/merging.h
    include/globfit2/optimization/mergingFunctors.h
    include/optimization/qp/gurobiOpt.h
    include/optimization/qp/MyGRBCallback.h
    include/optimization/qp/solver.h
)

SET( GLOBOPT_SRC_LIST
    src/qp_main.cpp
    src/subsample.cpp
    src/segment3D.cpp
    src/generate.cpp
    src/problemSetup.cpp
    src/solve.cpp
    src/merge.cpp
    src/datafit.cpp
    #src/optimization/candidateGenerator.cpp # cached compilation
    #src/optimization/qp/visualizer.cpp     # cached compilation
)
IF(WITH_SAMPLE_INPUT)
    SET( GLOBOPT_H_LIST ${GLOBOPT_H_LIST}
        include/globfit2/ground_truth/gtCreator.h
    )
    SET( GLOBOPT_SRC_LIST ${GLOBOPT_SRC_LIST}
        src/ground_truth/ground_truth.cpp
        src/ground_truth/gtCreator.cpp
    )
ENDIF(WITH_SAMPLE_INPUT)

ADD_EXECUTABLE( ${GLOBOPT_TARGET_NAME} ${GLOBOPT_SRC_LIST} ${GLOBOPT_HPP_LIST} ${GLOBOPT_H_LIST})
TARGET_LINK_LIBRARIES( ${GLOBOPT_TARGET_NAME}
    ${GUROBI_LIBRARIES}
    ${MOSEK_LIBRARIES}
    ${BONMIN_LIBRARIES}
#    ${GLOBOPTVIS_LIBRARIES}
#    ${WORKSPACE_DIR}/globOpt/globOpt/../visualization/build/Release/lib/libglobOptVis.a
    ${OpenCV_LIBRARIES}
    ${PCL_LIBRARIES}
    #pthread
    boost_filesystem
    boost_system
    boost_thread
)

message( "${PATH_VIS_ROOT}/build/Release/lib/libglobOptVis.a" )

## _______________________________________________________________________________________________
#SET( TEST_TARGET test_bst)
#add_executable( ${TEST_TARGET} src/test/test_main.cpp )
#TARGET_LINK_LIBRARIES( ${TEST_TARGET}
#    ${PCL_LIBRARIES}
#)
