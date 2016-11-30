cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(psi4pasture LANGUAGES C CXX VERSION 0.0.1)

# sortout transqt2/ccsort interdependence
if(ENABLE_transqt2 OR ENABLE_ccsort)
  if(NOT (ENABLE_ccsort AND ENABLE_transqt2))
    message(STATUS  "ccsort and transqt2 enabled since interdependent")
    set(ENABLE_ccsort "ON")
    set(ENABLE_transqt2 "ON")
  endif()
endif()
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
find_package(psi4 1.0 REQUIRED)

# define the install path
get_filename_component(psi4_lib_loc ${core_LIBRARY} DIRECTORY)
set(install_loc "${psi4_lib_loc}/driver/pasture")

message(STATUS "Will install pasture modules to ${psi4_lib_loc}/driver/pasture")
foreach(pasture_dir ccsort transqt2)
  if(${ENABLE_${pasture_dir}})
    add_subdirectory(${pasture_dir})
  endif()

endforeach()
install(FILES __init__.py addins.py locations.py
        DESTINATION ${install_loc})
