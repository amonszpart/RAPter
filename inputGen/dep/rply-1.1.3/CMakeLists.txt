set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(rply_ROOT "${dep_ROOT}/rply-1.1.3")

#############################################
set(rply_SRC
    ${rply_ROOT}/rply.c)
set(rply_INCLUDE
    ${rply_ROOT}/rply.h)

add_library(rply ${rply_SRC} ${rply_INCLUDE})
