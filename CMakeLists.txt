cmake_minimum_required(VERSION 3.16)
project(MARSQ LANGUAGES NONE)

# Preserve the registered native regressions while capturing the full source tree.
enable_testing()
add_subdirectory(tests)
