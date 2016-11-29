cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(psi4pasture LANGUAGES C CXX VERSION 0.0.1)

# sortout transqt2/ccsort interdependence
if(ENABLE_transqt2)
  if(NOT ENABLE_ccsort)
    message(WARNING "transqt2 enabled without ccsort: ccsort is required to use transqt2")
    message(STATUS  "Enabling ccsort")
    set(ENABLE_ccsort ON)
  endif()
elseif(ENABLE_ccsort)
  if(NOT ENABLE_transqt2)
    message(WARNING "ccsort enabled without transqt2: transqt2 is required to use ccsort")
    message(STATUS  "Enabling transqt2")
    set(ENABLE_transqt2 ON)
  endif()
endif()

find_package(psi4 1.0 REQUIRED)
# find the top level
get_filename_component(psi4_lib_loc ${core_LIBRARY} DIRECTORY)


set(install_loc "${psi4_lib_loc}/driver/pasture")
message(STATUS "Will install pasture modules to ${psi4_lib_loc}/driver/pasture")
foreach(pasture_dir ccsort transqt2)
  if(${ENABLE_${pasture_dir}})
    add_subdirectory(${pasture_dir})
  endif()
endforeach()
install(FILES __init__.py addins.py metadata.py
        DESTINATION ${install_loc})
