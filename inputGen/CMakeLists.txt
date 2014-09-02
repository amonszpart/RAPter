PROJECT(inputGen)
cmake_minimum_required(VERSION 2.8)

FIND_PACKAGE(Qt4 COMPONENTS QtCore QtGui QtOpenGL REQUIRED)

FIND_PACKAGE(OpenGL REQUIRED)
FIND_PACKAGE(Eigen3 REQUIRED)

SET(inputGen_SOURCES
  src/main.cpp
  src/mainwindow.cpp
  src/samplegeneratorfactory.cpp
  src/myview.cpp
  src/myscene.cpp)
SET(inputGen_HEADERS
  include/mainwindow.h
  include/samplegenerator.h
  include/samplegeneratorfactory.h
  include/myview.h
  include/myscene.h
  include/types.h
  include/primitive.h)
SET(inputGen_FORMS
  ui/mainwindow.ui
  ui/samplegeneratorfactory.ui)
# SET(inputGen_RESOURCES images.qrc)

SET(QT_USE_QTOPENGL TRUE)

QT4_WRAP_CPP(inputGen_HEADERS_MOC ${inputGen_HEADERS})
QT4_WRAP_UI(inputGen_FORMS_HEADERS ${inputGen_FORMS})

INCLUDE(${QT_USE_FILE})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include)
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${EIGEN3_INCLUDE_DIR}/unsupported)

ADD_DEFINITIONS(${QT_DEFINITIONS})

ADD_EXECUTABLE(inputGen ${inputGen_SOURCES}
    ${inputGen_HEADERS}
    ${inputGen_HEADERS_MOC}
    ${inputGen_FORMS_HEADERS} )
TARGET_LINK_LIBRARIES(inputGen ${OPENGL_LIBRARIES} ${QT_LIBRARIES})
