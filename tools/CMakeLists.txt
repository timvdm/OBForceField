set(QT_USE_QTOPENGL true)

include(${QT_USE_FILE})

# Ensure the Avogadro include directory is always first
get_directory_property(tmp_include_dirs INCLUDE_DIRECTORIES)
set_directory_properties(PROPERTIES INCLUDE_DIRECTORIES
  "${libavogadro_BINARY_DIR}/include;${tmp_include_dirs}")
include_directories(
  ${CMAKE_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${EIGEN2_INCLUDE_DIR}
  ${OPENBABEL2_INCLUDE_DIR}
)

link_directories(
  ${CMAKE_BINARY_DIR}
  ${OPENBABEL2_LIBRARY_DIRS}
)

set(tools
    energy
    minimize
)

foreach (tool ${tools})
  message(STATUS "Tool:  ${tool}")
  set(tool_SRCS ${tool}.cpp)
  add_executable(${tool} ${tool_SRCS})
  target_link_libraries(${tool}
    ${OPENBABEL2_LIBRARIES}
    ${QT_LIBRARIES}
    ${QT_QTTEST_LIBRARY}
    obforcefields)
endforeach (tool ${tools})
