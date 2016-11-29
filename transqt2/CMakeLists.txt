set(transqt2_sources transtwo_rhf.cc cache.cc get_params.cc transtwo_uhf.cc transqt.cc transone.cc get_moinfo.cc file_build_presort.cc idx_error.cc idx_permute_presort.cc )
add_library(transqt2 MODULE ${transqt2_sources})
target_include_directories(transqt2 SYSTEM PRIVATE ${PSI4_INCLUDE_DIRS})
target_link_libraries(transqt2 ${core_LIBRARY})
set_target_properties(transqt2 PROPERTIES
  PREFIX ""
  COMPILE_FLAGS "-std=c++11")
if(APPLE)
  set_target_properties(transqt2 PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif()

install(TARGETS transqt2
        LIBRARY DESTINATION ${install_loc})
