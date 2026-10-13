cmake_minimum_required(VERSION 3.16)
project(marsq_pamera_regression NONE)
enable_testing()
find_package(Python3 REQUIRED COMPONENTS Interpreter)
add_test(NAME pamera_singular_block
  COMMAND ${Python3_EXECUTABLE} -m unittest discover -s "${CMAKE_CURRENT_SOURCE_DIR}" -p "test_*.py" -v)
set_tests_properties(pamera_singular_block PROPERTIES TIMEOUT 120
  ENVIRONMENT "TMPDIR=${CMAKE_CURRENT_BINARY_DIR}")
