PROJECT(inputGen)
cmake_minimum_required(VERSION 2.8)

FIND_PACKAGE(Qt4 COMPONENTS QtCore QtGui QtXml QtOpenGL REQUIRED)

FIND_PACKAGE(OpenGL REQUIRED)
FIND_PACKAGE(Eigen3 REQUIRED)
#FIND_PACKAGE(Boost COMPONENTS random REQUIRED)

SET(inputGen_SOURCES
  src/main.cpp
  src/mainwindow.cpp
  src/mergedialog.cpp
  src/samplerfactory.cpp
  src/displacementfactory.cpp
  src/myview.cpp
  src/myscene.cpp)
SET(inputGen_HEADERS
  include/mainwindow.h
  include/mergedialog.h
  include/sampler.h
  include/samplerfactory.h
  include/displacement.h
  include/displacementfactory.h
  include/myview.h
  include/myscene.h
  include/types.h
  include/typesGL.h
  include/project.h
  include/primitive.h
  include/convexHull2D.h)
SET(inputGen_IMPL
  include/impl/randomdisplacement.hpp
  include/impl/biasdisplacement.hpp
  include/impl/sampler.hpp
  include/impl/convexHull2D.hpp)
SET(inputGen_FORMS
  ui/mainwindow.ui
  ui/mergedialog.ui
  ui/samplerfactory.ui
  ui/displacementfactory.ui)
# SET(inputGen_RESOURCES images.qrc)

SET(QT_USE_QTOPENGL TRUE)

QT4_WRAP_CPP(inputGen_HEADERS_MOC ${inputGen_HEADERS})
QT4_WRAP_UI(inputGen_FORMS_HEADERS ${inputGen_FORMS})

ADD_DEFINITIONS( -std=c++11 )

INCLUDE(${QT_USE_FILE})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include)
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/dep)
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR}/unsupported)

add_subdirectory("${PROJECT_SOURCE_DIR}/dep")

ADD_DEFINITIONS(${QT_DEFINITIONS})

ADD_EXECUTABLE(inputGen ${inputGen_SOURCES}
    ${inputGen_HEADERS}
    ${inputGen_IMPL}
    ${inputGen_HEADERS_MOC}
    ${inputGen_FORMS_HEADERS} )
TARGET_LINK_LIBRARIES(inputGen
    ${OPENGL_LIBRARIES}
    ${QT_LIBRARIES}
#    ${Boost_LIBRARIES}
    rply
)
