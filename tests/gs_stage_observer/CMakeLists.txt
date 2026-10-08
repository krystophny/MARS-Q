cmake_minimum_required(VERSION 3.20)
project(mars_gs_stage_observer LANGUAGES Fortran)
enable_testing()
set(CHEASE_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/../../CheaseMerge/chease.f" CACHE FILEPATH "Native observer source")
file(READ "${CHEASE_SOURCE}" NATIVE)
string(FIND "${NATIVE}" "      MODULE folded_band" END)
if(END LESS 0)
  message(FATAL_ERROR "Missing native observer module boundary")
endif()
string(SUBSTRING "${NATIVE}" 0 ${END} OBSERVER)
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/native_observer.f" "${OBSERVER}")
add_executable(observer_oracle "${CMAKE_CURRENT_BINARY_DIR}/native_observer.f" oracle.f90)
set_target_properties(observer_oracle PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
target_include_directories(observer_oracle PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND CMAKE_BUILD_TYPE STREQUAL "Debug")
  target_compile_options(observer_oracle PRIVATE -fcheck=all -ffpe-trap=invalid,zero,overflow)
endif()
add_test(NAME observer_native_headers_and_unchanged_arrays
  COMMAND ${CMAKE_COMMAND} -DPROGRAM=$<TARGET_FILE:observer_oracle>
    -DOUTDIR=${CMAKE_CURRENT_BINARY_DIR}/oracle-output -P ${CMAKE_CURRENT_SOURCE_DIR}/run.cmake)
add_test(NAME observer_disabled_no_array_mutation
  COMMAND ${CMAKE_COMMAND} -E env CHEASE_GS_OBSERVER_DIR= $<TARGET_FILE:observer_oracle> disabled)

string(FIND "${NATIVE}" "      subroutine GS_CAPTURE_MATRIX" BEGIN)
string(SUBSTRING "${NATIVE}" ${BEGIN} -1 TAIL)
string(FIND "${TAIL}" "\n      end\n" END)
math(EXPR END "${END}+11")
string(SUBSTRING "${TAIL}" 0 ${END} MATRIX_HOOK)
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/native_matrix_hook.f" "${MATRIX_HOOK}")
add_executable(hook_oracle "${CMAKE_CURRENT_BINARY_DIR}/native_observer.f"
  "${CMAKE_CURRENT_BINARY_DIR}/native_matrix_hook.f" hook_oracle.f)
set_target_properties(hook_oracle PROPERTIES Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/hook-modules")
target_include_directories(hook_oracle PRIVATE "${CMAKE_CURRENT_BINARY_DIR}/hook-modules"
  "${CMAKE_CURRENT_SOURCE_DIR}/../../CheaseMerge")
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  target_compile_options(hook_oracle PRIVATE -fdefault-real-8 -fdefault-double-8 -mcmodel=medium)
  target_link_options(hook_oracle PRIVATE -mcmodel=medium)
endif()
add_test(NAME observer_native_hook_sigma_not_profile_radius
  COMMAND ${CMAKE_COMMAND} -DPROGRAM=$<TARGET_FILE:hook_oracle>
    -DOUTDIR=${CMAKE_CURRENT_BINARY_DIR}/hook-output -P ${CMAKE_CURRENT_SOURCE_DIR}/run.cmake)
