cmake_minimum_required(VERSION 3.16)
project(marsq_regression NONE)
enable_testing()
find_package(Python3 REQUIRED COMPONENTS Interpreter)
set(MARS_NATIVE_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../MarsQ_2FK" CACHE PATH
    "Native MARS source inputs for compiled regression tests")
set(MARS_DIAGNOSTIC_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../tools" CACHE PATH
    "Passive MARS diagnostic tools used by regression tests")
set(MARS_PAMS_FORTSPARSE_PREFIX "$ENV{MARS_PAMS_FORTSPARSE_PREFIX}" CACHE PATH
    "Existing native fortsparse build for the optional PAMS/CSC regression")
set(MARS_PAMS_LINK_LIBRARIES "$ENV{MARS_PAMS_LINK_LIBRARIES}" CACHE STRING
    "Space-separated extra link library names for the PAMS/CSC regression")
file(GLOB tests CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/test_*.py")
foreach(t IN LISTS tests)
  get_filename_component(name "${t}" NAME_WE)
  if(name STREQUAL "test_pams_global_csc")
    set(pams_arguments --ctest --build-parent "${CMAKE_CURRENT_BINARY_DIR}")
    if(NOT "${MARS_PAMS_FORTSPARSE_PREFIX}" STREQUAL "")
      list(APPEND pams_arguments --fortsparse-prefix "${MARS_PAMS_FORTSPARSE_PREFIX}")
    endif()
    separate_arguments(pams_libraries NATIVE_COMMAND "${MARS_PAMS_LINK_LIBRARIES}")
    foreach(library IN LISTS pams_libraries)
      list(APPEND pams_arguments --link-library "${library}")
    endforeach()
    add_test(NAME ${name} COMMAND ${Python3_EXECUTABLE} "${t}" ${pams_arguments}
             WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
    set_tests_properties(${name} PROPERTIES SKIP_RETURN_CODE 77)
  else()
    add_test(NAME ${name} COMMAND ${Python3_EXECUTABLE} -m unittest -v ${name}
             WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}")
  endif()
  set_tests_properties(${name} PROPERTIES TIMEOUT 900 ENVIRONMENT "TMPDIR=${CMAKE_CURRENT_BINARY_DIR}")
endforeach()

# The snapshot regression compiles these production inputs in a private build.
# Declare its sibling source tree to CTest and immutable execution providers.
set(snapshot_inputs)
foreach(name anisotropic.f kinetic.f globalm.f dimensim.f rcomdm.f kineticm.f
             anisotropicm.f compam.inc)
  list(APPEND snapshot_inputs "${MARS_NATIVE_SOURCE_DIR}/${name}")
endforeach()
set_property(TEST test_energy_snapshot APPEND PROPERTY ENVIRONMENT
             "MARS_NATIVE_SOURCE_DIR=${MARS_NATIVE_SOURCE_DIR}")
set_property(TEST test_energy_snapshot PROPERTY REQUIRED_FILES "${snapshot_inputs}")

# Keep the offline oracle inputs explicit in frozen CTest execution views.
foreach(name test_energy_oracle test_energy_snapshot_check)
  if(TEST ${name})
    set_property(TEST ${name} APPEND PROPERTY ENVIRONMENT
                 "MARS_DIAGNOSTIC_SOURCE_DIR=${MARS_DIAGNOSTIC_SOURCE_DIR}")
    set_property(TEST ${name} PROPERTY REQUIRED_FILES
                 "${MARS_DIAGNOSTIC_SOURCE_DIR}/energy_oracle.py;${MARS_DIAGNOSTIC_SOURCE_DIR}/energy_snapshot_check.py")
  endif()
endforeach()
