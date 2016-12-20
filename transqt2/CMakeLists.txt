set(transqt2_sources transtwo_rhf.cc cache.cc get_params.cc transtwo_uhf.cc transqt.cc transone.cc get_moinfo.cc file_build_presort.cc idx_error.cc idx_permute_presort.cc )
add_psi4_plugin(transqt2 ${transqt2_sources})
target_include_directories(transqt2 SYSTEM PRIVATE ${psi4_INCLUDE_DIR})

install(TARGETS transqt2
        LIBRARY DESTINATION ${install_loc})
