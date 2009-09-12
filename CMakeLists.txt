project(OBForceFields) 
cmake_minimum_required(VERSION 2.6.0)

message(STATUS "OBForceFields 0.9.0")
message(STATUS "source: ${OBForceFields_SOURCE_DIR}")

set(CMAKE_MODULE_PATH ${OBForceFields_SOURCE_DIR}/cmake/modules)
find_package(OpenBabel2 REQUIRED) # find and setup OpenBabel
find_package(Eigen2 REQUIRED) # find and setup Eigen2

find_package(OpenCL REQUIRED)
link_directories(${OPENCL_LIBRARY_DIR})


# QtCore for QtConcurrent
find_package(Qt4) # find and setup Qt4 for this project
include(${QT_USE_FILE})

enable_testing()

include_directories(
    ${OBForceFields_SOURCE_DIR}/include
    ${OPENBABEL2_INCLUDE_DIR}
    ${EIGEN2_INCLUDE_DIR}
    ${QT_INCLUDES}
    ${CMAKE_BINARY_DIR}
    ${OPENCL_INCLUDE_DIR}
)

set(obforcefields_srcs
    src/obvariant.cpp
    src/obforcefield.cpp
    src/obfunction.cpp
    src/obfunctionterm.cpp
    src/obminimize.cpp
    src/oblogfile.cpp
    src/vectormath.cpp

    src/forcefields/mmff94/common.cpp
    src/forcefields/mmff94/parameter.cpp
    src/forcefields/mmff94/function.cpp
    src/forcefields/mmff94/bond.cpp
    src/forcefields/mmff94/angle.cpp
    src/forcefields/mmff94/strbnd.cpp
    src/forcefields/mmff94/torsion.cpp
    src/forcefields/mmff94/oop.cpp
    src/forcefields/mmff94/vdw.cpp
    src/forcefields/mmff94/electro.cpp
    
    src/forcefields/mmff94/vdw_opencl.cpp
    src/forcefields/mmff94/electro_opencl.cpp
)

add_library(obforcefields SHARED ${obforcefields_srcs})
target_link_libraries(obforcefields 
    ${OPENBABEL2_LIBRARIES}
    ${OPENCL_LIBRARIES}
)


add_subdirectory(tests)
add_subdirectory(tools)

########################################
#  Qt4
########################################
#if (QT_FOUND)
  #message(STATUS "Found Qt4 libaries -- Enabling QtPainter.")
  #set(qtdepict_srcs ${depict_srcs} src/qtpainter.cpp 
      #src/qtdepict/qtdepict.cpp src/qtdepict/mainwindow.cpp)
#
  #qt4_wrap_ui(mainwindow_UIS src/qtdepict/mainwindow.ui)
  #qt4_wrap_cpp(mainwindow_MOC src/qtdepict/mainwindow.h)

#  add_executable(qtdepict ${mainwindow_UIS} ${qtdepict_srcs} ${mainwindow_MOC})
  #target_link_libraries(qtdepict ${OPENBABEL2_LIBRARIES} ${QT_LIBRARIES})
  #install(TARGETS qtdepict DESTINATION bin)
#endif (QT_FOUND)


