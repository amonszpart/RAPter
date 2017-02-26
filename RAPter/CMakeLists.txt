SET(PROJECT_NAME RAPter)
PROJECT(${PROJECT_NAME})
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#_____________________________________#
#_______________PATHS_________________#
#_____________________________________#

SET( WORKSPACE_DIR "$ENV{HOME}/workspace" )
SET( THIRD_PARTY_DIR "${WORKSPACE_DIR}/3rdparty" )

# PCL
SET( PCL_DIR "${THIRD_PARTY_DIR}/pcl-trunk2/install/share/pcl-1.8/"
        CACHE FILEPATH "Folder containing \"PCLConfig.cmake\", usually \"{PCL_INSTALL_DIR}/share/pcl-x.x/\"." )
# BONMIN
SET( PATH_BONMIN_DIR "${THIRD_PARTY_DIR}/CoinBonmin-stable/build/"
    CACHE FILEPATH "Folder containing the compiled Bonmin libraries")
#SET( PATH_BONMIN_DIR "${THIRD_PARTY_DIR}/Bonmin-1.8.0-lapack/")

## QCQPCPP
#SET( PATH_QCQPCPP_ROOT "${WORKSPACE_DIR}/QCQPcpp/"
#        CACHE FILEPATH "Root directory of QCQPcpp library containing the \"include\" directory.")

# PEARL, GCO
#SET( PATH_GLOBOPT_DEPS_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../../globOptDeps"
#        CACHE PATH "Folder containing additional dependency libraries (pearl, gco-v3, etc.)" FORCE )

# FBI
#SET( PATH_LIBFBI_ROOT "${THIRD_PARTY_DIR}/libfbi/")

#_____________________________________#
#______________TARGETS________________#
#_____________________________________#

SET( WITH_PEARL OFF CACHE BINARY "Compile pearl executable.")
SET( WITH_RSAC OFF CACHE BINARY "Compile ransac executable.")
SET( WITH_CORRESP OFF CACHE BINARY "Correspondance miniprogram compilation.")
SET( WITH_GLOBFIT OFF CACHE BINARY "GlobFit (Li'11) comparison.")
SET( WITH_EVAL OFF CACHE BINARY "Evaluation code")
SET( WITH_REFIT OFF CACHE BINARY "Constrained refit project")
SET( WITH_GCO OFF CACHE BINARY "Compile alpha-expansion library by Veksler and Delong, needed by PEARL, RSAC and REFIT projects.")
#SET( WITH_GAUSSSPHERE OFF CACHE BINARY "Compile gaussSphere." )
SET( WITH_TO_PS OFF CACHE BINARY "Compile primitives to ps converter." )
#SET( WITH_PLYCONVERTER ON CACHE BINARY "Compile ply-converter executable.")

#_____________________________________#
#________________BUILD________________#
#_____________________________________#

SET( dir ${CMAKE_CURRENT_SOURCE_DIR}/build/${CMAKE_BUILD_TYPE} )
SET( EXECUTABLE_OUTPUT_PATH ${dir}/bin CACHE PATH "Build directory" FORCE)
SET( LIBRARY_OUTPUT_PATH ${dir}/lib CACHE PATH "Build directory" FORCE)

SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -Wall -Wno-unused-local-typedefs " ) #-Wreturn-type
SET( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE} -O2" )
ADD_DEFINITIONS(-DRAPTER_WITH_BONMIN)

#____________________________________#
#________________DEPS________________#
#____________________________________#

# OpenGL
FIND_PACKAGE( OpenGL REQUIRED )
#INCLUDE_DIRECTORIES( ${OPENGL_INCLUDE_DIRS} )

# PCL
FIND_PACKAGE( PCL REQUIRED ) # COMPONENTS COMMON IO SEARCH FEATURES FILTERS SEGMENTATION VISUALIZATION SAMPLE_CONSENSUS
#INCLUDE_DIRECTORIES( ${PCL_INCLUDE_DIRS} )
#ADD_DEFINITIONS( ${PCL_DEFINITIONS} )

# OpenMP
FIND_PACKAGE(OpenMP)
IF(OPENMP_FOUND)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
ENDIF(OPENMP_FOUND)

# Bonmin
IF(EXISTS ${PATH_BONMIN_DIR}/lib/libcoinhsl.so)
  SET (BONMIN_SOLVER_LIB ${PATH_BONMIN_DIR}/lib/libcoinhsl.so)
  MESSAGE(STATUS "HSL Solver available, using HSL instead of MUMPS")
ELSE()
    SET (BONMIN_SOLVER_LIB ${PATH_BONMIN_DIR}/lib/libcoinmumps.so)
ENDIF(EXISTS ${PATH_BONMIN_DIR}/lib/libcoinhsl.so)

SET(BONMIN_LIBRARIES
    ${PATH_BONMIN_DIR}/lib/libipopt.so
    ${PATH_BONMIN_DIR}/lib/libbonmin.so
    ${PATH_BONMIN_DIR}/lib/libCbc.so
    ${PATH_BONMIN_DIR}/lib/libCgl.so
    ${PATH_BONMIN_DIR}/lib/libClp.so
    ${BONMIN_SOLVER_LIB}
    ${PATH_BONMIN_DIR}/lib/libCoinUtils.so
    ${PATH_BONMIN_DIR}/lib/libOsi.so
    ${PATH_BONMIN_DIR}/lib/libOsiClp.so
    lapack
    dl
)
SET( BONMIN_INCLUDE_DIRS "${PATH_BONMIN_DIR}/include" )

# QCQPCPP
SET( QCQPCPP_H_LIST
    external/qcqpcpp/include/qcqpcpp/optProblem.h
    external/qcqpcpp/include/qcqpcpp/bonminOptProblem.h
    external/qcqpcpp/include/qcqpcpp/io/io.h
)
SET( QCQPCPP_HPP_LIST
    external/qcqpcpp/include/qcqpcpp/impl/optProblem.hpp
)
SET( QCQPCPP_INCLUDE_DIRS "external/qcqpcpp/include" )

# AlphaExpansion (GCO lib)
IF(WITH_GCO OR WITH_PEARL OR WITH_REFIT)
    add_definitions( -DGCO_ENERGYTYPE=float )
    add_subdirectory( external/gco-v3-float )
    #SET( PATH_GCO_ROOT ${PATH_GLOBOPT_DEPS_ROOT}/gco-v3-float CACHE PATH "alpha-expansion root" FORCE )

    #add_library( libgco-v3 SHARED IMPORTED )
    #set_property( TARGET libgco-v3 PROPERTY IMPORTED_LOCATION_RELEASE ${LIBRARY_OUTPUT_PATH}/libgco-v3-0-float.so )
    #set_property( TARGET libgco-v3 PROPERTY IMPORTED_LOCATION_DEBUG   ${LIBRARY_OUTPUT_PATH}/libgco-v3-0-floatd.so )
    SET( GCO_INCLUDE_DIRS external/gco-v3-float/include )
    #SET( GCO_LIBRARIES "libgco-v3" )
    SET( GCO_LIBRARIES "gco-v3-0-float" )
ENDIF()

# libfbi
#INCLUDE_DIRECTORIES( ${PATH_LIBFBI_ROOT}/include )

#_______________________________________#
#____________TemplateInst_______________#
#_______________________________________#

# Common template instantiations for all projects
SET( TEMPLATE_INST_SRC_LIST
    src/templateInstantiation/taggable.cpp
    src/templateInstantiation/primitive.cpp
    src/templateInstantiation/linePrimitive.cpp
    src/templateInstantiation/planePrimitive.cpp
    src/templateInstantiation/segmentation.cpp
    src/templateInstantiation/assignmentOps.cpp
)

##___________________________________________________________________________##
##                                  RAPTER                                   ##
##___________________________________________________________________________##
SET(RAPTER_TARGET_NAME rapter)

ADD_DEFINITIONS( -DRAPTER_USE_PCL=1
                 -DVTK_EXCLUDE_STRSTREAM_HEADERS
               )

SET( RAPTER_HPP_LIST
    include/rapter/io/impl/io.hpp
    include/rapter/io/inputParser.hpp
    include/rapter/optimization/impl/segmentation.hpp
    include/rapter/optimization/impl/solver.hpp
    include/rapter/optimization/impl/problemSetup.hpp
    include/rapter/optimization/impl/merging.hpp
    include/rapter/optimization/impl/candidateGenerator.hpp
    include/rapter/primitives/impl/taggable.hpp
    include/rapter/primitives/impl/planePrimitive.hpp
    include/rapter/primitives/impl/linePrimitive.hpp
    include/rapter/processing/util.hpp
    include/rapter/processing/impl/angleUtil.hpp
    include/rapter/processing/graph.hpp
    include/rapter/processing/diagnostic.hpp
    include/rapter/processing/impl/angle.hpp
    include/rapter/util/diskUtil.hpp
    include/rapter/util/util.hpp
    include/rapter/util/impl/pclUtil.hpp
    include/rapter/util/containers.hpp
    include/rapter/util/lruCache.hpp
    ${QCQPCPP_HPP_LIST}
)

SET( RAPTER_H_LIST
    include/rapter/simpleTypes.h
    include/rapter/typedefs.h
    include/rapter/parameters.h
    include/rapter/io/io.h
    include/rapter/optimization/patchDistanceFunctors.h
    include/rapter/optimization/energyFunctors.h
    include/rapter/optimization/candidateGenerator.h
    include/rapter/optimization/problemSetup.h
    include/rapter/optimization/merging.h
    include/rapter/optimization/mergingFunctors.h
    include/rapter/optimization/segmentation.h
    include/rapter/optimization/solver.h
    include/rapter/primitives/angles.h
    include/rapter/primitives/linePrimitive.h
    include/rapter/primitives/taggable.h
    include/rapter/primitives/primitive.h
    include/rapter/primitives/pointPrimitive.h
    include/rapter/primitives/planePrimitive.h
    include/rapter/util/parse.h
    include/rapter/util/pclUtil.h
    ${QCQPCPP_H_LIST}
)

SET( RAPTER_SRC_LIST
    src/main.cpp
    src/subsample.cpp
    src/segment.cpp
    src/generate.cpp
    src/generate3D.cpp
    src/problemSetup.cpp
    src/problemSetup3D.cpp
    src/solve.cpp
    src/solve3D.cpp
    src/merge.cpp
#    src/datafit.cpp
#    src/reassign.cpp
    src/represent.cpp
    ${TEMPLATE_INST_SRC_LIST}
)

INCLUDE_DIRECTORIES(
    include
    ${QCQPCPP_INCLUDE_DIRS}
    ${BONMIN_INCLUDE_DIRS}
    ${PCL_INCLUDE_DIRS}
    ${EIGEN_INCLUDE_DIRS}
)

ADD_DEFINITIONS( ${PCL_DEFINITIONS} )

ADD_EXECUTABLE( ${RAPTER_TARGET_NAME}
    ${RAPTER_H_LIST}
    ${RAPTER_HPP_LIST}
    ${RAPTER_SRC_LIST}
)

TARGET_LINK_LIBRARIES( ${RAPTER_TARGET_NAME}
    ${BONMIN_LIBRARIES}
    ${PCL_LIBRARIES}
    boost_filesystem
    boost_system
    boost_thread
)

##___________________________________________________________________________##
##                                  PEaRL                                    ##
##___________________________________________________________________________##

IF(WITH_PEARL)
    SET( PEARL_TARGET pearl)

#    find_package(OpenMP)
#    if (OPENMP_FOUND)
#        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    endif()

    ## Deps
    #SET( PATH_PEARL_ROOT ${PATH_GLOBOPT_DEPS_ROOT}/pearl CACHE PATH "Pearl root" FORCE )
    #SET( PATH_PEARL_ROOT pearl)
    #SET( WITH_GCO ON CACHE BINARY "Use alpha-expansion library" FORCE)
    #SET( PATH_GCO_ROOT ${PATH_GLOBOPT_DEPS_ROOT}/gco-v3-float CACHE PATH "alpha-expansion root" FORCE )
    IF( NOT ${WITH_GCO} )
        message( FATAL_ERROR "\nPlease enable WITH_GCO flag if you want to use the pearl project!\n\n" )
    ENDIF()

    ## Sources
    SET( PEARL_TARGET_H
        include/pearl/pearl.h
        #${PATH_GLOBOPT_DEPS_ROOT}/pearl/include/pearl/pearl.h
    )

    SET( PEARL_TARGET_HPP
        include/pearl/pearl.hpp
        include/rapter/comparison/impl/reassign.hpp
    )

    SET( PEARL_TARGET_SRC
        src/comparison/pearl.cpp
        src/reassign.cpp
        ${TEMPLATE_INST_SRC_LIST}
    )

    ## Linking
    INCLUDE_DIRECTORIES(
        ${GCO_INCLUDE_DIRS}
        ${PCL_INCLUDE_DIRS}
        ${EIGEN_INCLUDE_DIRS}
    )

    ADD_DEFINITIONS( ${PCL_DEFINITIONS} )

    ADD_EXECUTABLE( ${PEARL_TARGET}
        ${PEARL_TARGET_SRC}
        ${PEARL_TARGET_H}
        ${PEARL_TARGET_HPP}
    )

    TARGET_LINK_LIBRARIES( ${PEARL_TARGET}
        ${GCO_LIBRARIES}
        ${PCL_LIBRARIES}
        boost_filesystem
        boost_system
        boost_thread
    )
ENDIF(WITH_PEARL)

##___________________________________________________________________________##
##                                  RSAC                                     ##
##___________________________________________________________________________##

IF(WITH_RSAC)
    SET( RSAC_TARGET ransac)
    SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unused-local-typedefs " ) #-Wreturn-type # removed c++11

    #################################
    # libSchnabel
    set( SCHNABEL07_DIR external/schnabel07 )
    add_subdirectory( ${SCHNABEL07_DIR} )
    #add_library( libSchnabel07 SHARED IMPORTED )
    #set_property( TARGET libSchnabel07 PROPERTY IMPORTED_LOCATION_RELEASE ${LIBRARY_OUTPUT_PATH}/libSchnabel07.so  )
    #set_property( TARGET libSchnabel07 PROPERTY IMPORTED_LOCATION_DEBUG   ${LIBRARY_OUTPUT_PATH}/libSchnabel07d.so )
    #add_library( libSchnabel07_MiscLib SHARED IMPORTED )
    #set_property( TARGET libSchnabel07_MiscLib PROPERTY IMPORTED_LOCATION_RELEASE ${LIBRARY_OUTPUT_PATH}/libMiscLib.so  )
    #set_property( TARGET libSchnabel07_MiscLib PROPERTY IMPORTED_LOCATION_DEBUG   ${LIBRARY_OUTPUT_PATH}/libMiscLibd.so )
    #MESSAGE(STATUS "LIBRARY_OUTPUT_PATH: ${LIBRARY_OUTPUT_PATH}" )
    SET( SCHNABEL07_INCLUDE_DIRS ${SCHNABEL07_DIR} )
    #SET( SCHNABEL07_LIBRARIES libSchnabel07 libSchnabel07_MiscLib )
    SET( SCHNABEL07_LIBRARIES Schnabel07 )
    #################################

    ## Sources
    SET( RSAC_H
    )

    SET( RSAC_HPP
    )

    SET( RSAC_SRC
        src/comparison/ransac.cpp
        src/schnabelEnv.cpp
        ${TEMPLATE_INST_SRC_LIST}
    )

    ## Linking
    INCLUDE_DIRECTORIES(
        ${SCHNABEL07_INCLUDE_DIRS}
        ${PCL_INCLUDE_DIRS}
    )

    ADD_DEFINITIONS( ${PCL_DEFINITIONS} )

    ADD_EXECUTABLE( ${RSAC_TARGET}
        ${RSAC_H}
        ${RSAC_HPP}
        ${RSAC_SRC}
    )

    TARGET_LINK_LIBRARIES( ${RSAC_TARGET}
        ${SCHNABEL07_LIBRARIES}
        ${GCO_LIBRARIES}
        ${PCL_LIBRARIES}
        boost_filesystem
        boost_system
        boost_thread
    )
ENDIF(WITH_RSAC)

##___________________________________________________________________________##
##                           CORRESPONDENCE                                  ##
##___________________________________________________________________________##

IF(WITH_CORRESP)
    SET( CORRESP_TARGET "corresp" )

    SET( CORRESP_SRC
        src/correspondence.cpp
        ${TEMPLATE_INST_SRC_LIST}
    )

    INCLUDE_DIRECTORIES(
        ${PCL_INCLUDE_DIRS}
    )

    ADD_EXECUTABLE( ${CORRESP_TARGET}
        ${CORRESP_SRC}
        ${CORRESP_H}
    )

    TARGET_LINK_LIBRARIES( ${CORRESP_TARGET}
        ${PCL_LIBRARIES}
        boost_filesystem
    )
ENDIF(WITH_CORRESP)

##___________________________________________________________________________##
##______________________________toGlobFit____________________________________##
##___________________________________________________________________________##

IF(WITH_GLOBFIT)
    SET( GLOBFIT_TARGET "toGlobFit" )

    SET( GLOBFIT_H
        include/rapter/processing/assignmentOps.h
        include/rapter/comparison/globFit.h
        include/rapter/processing/subsamplePrimitives.h
    )

    SET( GLOBFIT_HPP
        include/rapter/processing/impl/assignmentOps.hpp
        include/rapter/comparison/impl/globFit.hpp
        include/rapter/processing/impl/subsamplePrimitives.hpp
    )

    SET( GLOBFIT_SRC
        src/comparison/globFit.cpp
        src/comparison/globFitTemplateInst.cpp
        src/comparison/subsampleTemplateInst.cpp
        ${TEMPLATE_INST_SRC_LIST}
    )

    # Linking
    INCLUDE_DIRECTORIES(
        ${PCL_INCLUDE_DIRS}
    )

    ADD_DEFINITIONS( ${PCL_DEFINITIONS} )

    ADD_EXECUTABLE( ${GLOBFIT_TARGET}
        ${GLOBFIT_H}
        ${GLOBFIT_HPP}
        ${GLOBFIT_SRC}
    )

    TARGET_LINK_LIBRARIES( ${GLOBFIT_TARGET}
        ${PCL_LIBRARIES}
        boost_filesystem
    )
ENDIF(WITH_GLOBFIT)

##___________________________________________________________________________##
##________________________________Eval_______________________________________##
##___________________________________________________________________________##

IF(WITH_EVAL)
    SET( EVAL_TARGET "eval" )

    SET( EVAL_H
        include/rapter/evaluation/assignPointsToTriangles.h
        include/rapter/primitives/planePrimitive.h
        include/rapter/primitives/linePrimitive.h
        include/rapter/io/trianglesFromObj.h
    )

    SET( EVAL_HPP
        include/rapter/evaluation/impl/assignPointsToTriangles.hpp
        include/rapter/primitives/impl/triangle.hpp
        include/rapter/primitives/impl/linePrimitive.hpp
        include/rapter/primitives/impl/planePrimitive.hpp
        include/rapter/io/ps.hpp
    )

    SET( EVAL_SRC
        src/evaluation/eval.cpp
        src/evaluation/assignPointsToTriangles.cpp
        src/templateInstantiation/trianglesFromObj.cpp
        ${TEMPLATE_INST_SRC_LIST}
    )

    ADD_EXECUTABLE( ${EVAL_TARGET}
        ${EVAL_SRC}
        ${EVAL_H}
        ${EVAL_HPP}
    )
    TARGET_LINK_LIBRARIES( ${EVAL_TARGET}
        ${PCL_LIBRARIES}
        boost_filesystem
    )
ENDIF(WITH_EVAL)

##___________________________________________________________________________##
##________________________________ReFit______________________________________##
##___________________________________________________________________________##

IF(WITH_REFIT)
    SET( REFIT_TARGET "refit" )


    SET( REFIT_H
        include/rapter/primitives/planePrimitive.h
        include/rapter/typedefs.h
        include/rapter/util/parse.h
        include/rapter/io/io.h
    )

    SET( REFIT_HPP
        include/rapter/primitives/impl/planePrimitive.hpp
        include/rapter/primitives/impl/linePrimitive.hpp
        include/rapter/processing/impl/angleUtil.hpp
    )

    SET( REFIT_SRC
        src/refit.cpp
        src/templateInstantiation/planePrimitive.cpp
        src/templateInstantiation/linePrimitive.cpp
        src/templateInstantiation/taggable.cpp
    )

    INCLUDE_DIRECTORIES(
        ${QCQPCPP_INCLUDE_DIRS}
        ${BONMIN_INCLUDE_DIRS}
        ${GCO_INCLUDE_DIRS}
        ${PCL_INCLUDE_DIRS}
    )

    ADD_DEFINITIONS( ${PCL_DEFINITIONS} )

    ADD_EXECUTABLE( ${REFIT_TARGET}
        ${REFIT_SRC}
        ${REFIT_H}
        ${REFIT_HPP}
    )
    TARGET_LINK_LIBRARIES( ${REFIT_TARGET}
        ${BONMIN_LIBRARIES}
        ${GCO_LIBRARIES}
        ${PCL_LIBRARIES}
        boost_filesystem
    )
ENDIF(WITH_REFIT)

##___________________________________________________________________________##
##                              PLY CONVERTER                                ##
##___________________________________________________________________________##

#IF(WITH_PLYCONVERTER)
#    SET( PLYCONVERTER_TARGET "convert2ply" )
#    SET( PLYCONVERTER_SRC
#        src/comparison/plyConverter.cpp
#    )

#    add_executable( ${PLYCONVERTER_TARGET} ${PLYCONVERTER_SRC} ${PLYCONVERTER_H})
#    TARGET_LINK_LIBRARIES( ${PLYCONVERTER_TARGET}
#        opencv_calib3d
#        opencv_highgui
#        opencv_imgproc
#        opencv_core
#        ${PCL_LIBRARIES}
#        boost_filesystem
#    )
#ENDIF(WITH_PLYCONVERTER)

##___________________________________________________________________________##
##______________________________GAUSSSPHERE__________________________________##
##___________________________________________________________________________##

#IF(WITH_GAUSSSPHERE)
#    SET( GAUSSSPHERE_TARGET "gaussSphere" )
#    SET( GAUSSSPHERE_SRC
#        src/gaussSphere.cpp
#    )

#    add_executable( ${GAUSSSPHERE_TARGET} ${GAUSSSPHERE_SRC} ${GAUSSSPHERE_H})
#    TARGET_LINK_LIBRARIES( ${GAUSSSPHERE_TARGET}
#        opencv_calib3d
#        opencv_highgui
#        opencv_imgproc
#        opencv_core
#        ${PCL_LIBRARIES}
#        boost_filesystem
#    )
#ENDIF(WITH_GAUSSSPHERE)

##___________________________________________________________________________##
##__________________________________TO_PS____________________________________##
##___________________________________________________________________________##

IF(WITH_TO_PS)
    SET( TO_PS_TARGET "toPs" )

    SET( TO_PS_H
    )

    SET( TO_PS_HPP
        include/rapter/io/ps.hpp
    )

    SET( TO_PS_SRC
        src/drawGraphs.cpp
        src/templateInstantiation/taggable.cpp
        src/templateInstantiation/primitive.cpp
        src/templateInstantiation/linePrimitive.cpp
        src/templateInstantiation/planePrimitive.cpp
    )

    ADD_EXECUTABLE( ${TO_PS_TARGET}
        ${TO_PS_SRC}
        ${TO_PS_H}
    )

    TARGET_LINK_LIBRARIES( ${TO_PS_TARGET}
        ${PCL_LIBRARIES}
        boost_filesystem
    )
ENDIF(WITH_TO_PS)
