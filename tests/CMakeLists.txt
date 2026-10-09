cmake_minimum_required(VERSION 3.16)
project(chease_regression NONE)
enable_testing()
find_package(Python3 REQUIRED COMPONENTS Interpreter)
get_filename_component(DEFAULT_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/.." ABSOLUTE)
set(CHEASE_ROOT "${DEFAULT_ROOT}" CACHE PATH "CHEASE source tree under test")
set(CHEASE_EXECUTABLE "" CACHE FILEPATH "Native Make-built CHEASE executable")
set(CHEASE_REFERENCE_EXECUTABLE "" CACHE FILEPATH "Reference executable for exact output comparison")
add_test(NAME regression COMMAND ${Python3_EXECUTABLE} -m unittest discover
  -s "${CMAKE_CURRENT_SOURCE_DIR}" -p "test_*.py" -v)
set_tests_properties(regression PROPERTIES TIMEOUT 300 ENVIRONMENT
  "CHEASE_ROOT=${CHEASE_ROOT};CHEASE_EXECUTABLE=${CHEASE_EXECUTABLE};CHEASE_REFERENCE_EXECUTABLE=${CHEASE_REFERENCE_EXECUTABLE};TMPDIR=${CMAKE_CURRENT_BINARY_DIR};OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
