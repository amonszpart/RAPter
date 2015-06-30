PROJECT(gco-v3-0-float)
cmake_minimum_required(VERSION 2.8)

# BUILD
#SET( dir ${CMAKE_CURRENT_SOURCE_DIR}/build )
#SET( EXECUTABLE_OUTPUT_PATH ${dir}/bin CACHE PATH "Build directory" FORCE)
#SET( LIBRARY_OUTPUT_PATH ${dir}/lib CACHE PATH "Build directory" FORCE)

add_definitions( -DGCO_ENERGYTYPE=float )

SET( HPP_LIST
    ${CMAKE_CURRENT_SOURCE_DIR}/include/gco/impl/graph.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/include/gco/impl/maxflow.hpp
)

SET( H_LIST
    ${CMAKE_CURRENT_SOURCE_DIR}/include/gco/block.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/gco/energy.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/gco/graph.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/gco/LinkedBlockList.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/gco/GCoptimization.h
)

INCLUDE_DIRECTORIES(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

SET( SRC_LIST
    ${CMAKE_CURRENT_SOURCE_DIR}/src/example.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/GCoptimization.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/LinkedBlockList.cpp
)

#add_executable( ${PROJECT_NAME} ${SRC_LIST} ${H_LIST} )
ADD_LIBRARY( ${PROJECT_NAME} SHARED ${SRC_LIST} ${H_LIST} ${HPP_LIST} )
SET_TARGET_PROPERTIES( ${PROJECT_NAME} PROPERTIES DEBUG_POSTFIX "d")

add_executable( test${PROJECT_NAME} ${SRC_LIST} ${H_LIST} )
