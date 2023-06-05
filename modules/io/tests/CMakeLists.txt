set(OST_IO_UNIT_TESTS
  test_io_pdb.py
  test_io_mmcif.py
  test_io_omf.py
  test_io_sdf.py
  test_clustal.cc
  test_io_pdb.cc
  test_io_crd.cc
  test_io_dcd.cc
  test_io_sdf.cc
  test_io_sdf_v3000.cc
  test_io_sequence_profile.cc
  test_pir.cc
  test_iomanager.cc
  tests.cc
  test_star_parser.cc
  test_mmcif_reader.cc
  test_mmcif_info.cc
  test_io_img.cc
  test_exceptions.cc
)

ost_unittest(MODULE io
             SOURCES "${OST_IO_UNIT_TESTS}"
             LINK ost_mol ost_seq)

if(NOT ENABLE_STATIC)
  add_executable(test_mae_standalone test_mae_standalone.cc)
  target_link_libraries(test_mae_standalone ost_mol)
  target_link_libraries(test_mae_standalone ost_io)
  target_link_libraries(test_mae_standalone ${Boost_REGEX_LIBRARY})
endif()
