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

# Disabled hydrogenscommand test for now - needs some work. Moc issues here.
# The commands should not just be included, they either need linking to
# or building. As plugin code is not part of the library it may require a
# different testing strategy.
set(tests
#  forcefield
  variant
  gradient
  mmff94parameterdb
  mmff94function
)

foreach (test ${tests})
  message(STATUS "Test:  ${test}")
  set(test_SRCS ${test}test.cpp)
  add_executable(${test}test ${test_SRCS})
  target_link_libraries(${test}test
    ${OPENBABEL2_LIBRARIES}
    ${QT_LIBRARIES}
    ${QT_QTTEST_LIBRARY}
    obforcefields)
  add_test(${test}Test ${CMAKE_BINARY_DIR}/tests/${test}test)
endforeach (test ${tests})
