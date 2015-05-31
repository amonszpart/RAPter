cmake_minimum_required(VERSION 2.6)
project(MiscLib)

aux_source_directory( ${CMAKE_CURRENT_SOURCE_DIR} SRC_LIST )

add_library( ${PROJECT_NAME} SHARED ${SRC_LIST} )
SET_TARGET_PROPERTIES( ${PROJECT_NAME} PROPERTIES DEBUG_POSTFIX "d")




