project(Schnabel07)
cmake_minimum_required(VERSION 2.6)

SET( CMAKE_CXX_FLAGS "-std=c++0x -Wall" ) #-Wreturn-type
#SET(CMAKE_VERBOSE_MAKEFILE true)
add_definitions( -DDOPARALLEL )

#SET( dir ${CMAKE_CURRENT_SOURCE_DIR}/build/${CMAKE_BUILD_TYPE})
#SET( EXECUTABLE_OUTPUT_PATH ${dir} CACHE PATH "Build directory" FORCE)
#SET( LIBRARY_OUTPUT_PATH    ${dir} CACHE PATH "Build directory" FORCE)

#ADD_SUBDIRECTORY(GfxTL)
ADD_SUBDIRECTORY(MiscLib)

aux_source_directory( . SRC_LIST )

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
MESSAGE(STATUS "cxxflags: " ${CMAKE_CXX_FLAGS})

SET( H_LIST
    ${H_LIST}
    GfxTL/AABox.h
    GfxTL/AABox.hpp
    GfxTL/AACube.h
    GfxTL/AACube.hpp
    GfxTL/AACubeTree.h
    GfxTL/AACubeTree.hpp
    GfxTL/AAKdTree.h
    GfxTL/AAKdTree.hpp
    GfxTL/KdTree.h
    GfxTL/KdTree.hpp
    GfxTL/AAPlane.h
    GfxTL/AAPlane.hpp
    GfxTL/BaseTree.h
    GfxTL/BaseTree.hpp
    GfxTL/NullStrategy.h
    GfxTL/NullStrategy.hpp
    GfxTL/Plane.h
    GfxTL/Plane.hpp
)

INCLUDE_DIRECTORIES( . )

add_library( ${PROJECT_NAME} SHARED ${SRC_LIST} ${H_LIST} )
SET_TARGET_PROPERTIES( ${PROJECT_NAME} PROPERTIES DEBUG_POSTFIX "d")

target_link_libraries(${PROJECT_NAME}
    MiscLib
    gomp )

# LIB COPY
#get_target_property( OUTPUT_LOCATION ${PROJECT_NAME} LOCATION )
#add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
#                   COMMAND ${CMAKE_COMMAND} -E copy_directory
#                   ${LIBRARY_OUTPUT_PATH}
#                   "${CMAKE_CURRENT_BINARY_DIR}/../lib/"
#)
