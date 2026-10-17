cmake_minimum_required(VERSION 3.16)
project(marsq_regression NONE)
enable_testing()
find_package(Python3 REQUIRED COMPONENTS Interpreter)
set(MARS_FORTSPARSE_PREFIX "" CACHE PATH
    "Existing fortsparse installation for the global sparse PAMS qualification")
if(MARS_FORTSPARSE_PREFIX)
  foreach(required IN ITEMS include lib/libfortsparse.a bin/fortsparse_umfpack_helper)
    if(NOT EXISTS "${MARS_FORTSPARSE_PREFIX}/${required}")
      message(FATAL_ERROR "MARS_FORTSPARSE_PREFIX is missing ${required}")
    endif()
  endforeach()
endif()
file(GLOB tests CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/test_*.py")
foreach(t IN LISTS tests)
  get_filename_component(name "${t}" NAME_WE)
  if(name STREQUAL "test_pams_global_csc")
    if(NOT MARS_FORTSPARSE_PREFIX)
      message(STATUS "Skipping test_pams_global_csc: set MARS_FORTSPARSE_PREFIX to qualify the optional global sparse backend")
      continue()
    endif()
    add_test(NAME ${name} COMMAND ${Python3_EXECUTABLE} "${t}"
             --fortsparse-prefix "${MARS_FORTSPARSE_PREFIX}"
             --build-parent "${CMAKE_CURRENT_BINARY_DIR}"
             --baseline-revision 8824bb18e1514b4a27f6357d5fe859d4fe690542
             --report "${CMAKE_CURRENT_BINARY_DIR}/global_caller_ctest_report.json")
  else()
    add_test(NAME ${name} COMMAND ${Python3_EXECUTABLE} -m unittest -v ${name}
             WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
  endif()
  set_tests_properties(${name} PROPERTIES TIMEOUT 900 ENVIRONMENT "TMPDIR=${CMAKE_CURRENT_BINARY_DIR}")
endforeach()
