set(OST_MOL_ALG_UNIT_TESTS
  test_superposition.cc
  tests.cc
  test_consistency_checks.cc
  test_partial_sec_struct_assignment.cc
  test_pdbize.py
  test_convenient_superpose.py
  test_hbond.py
  test_accessibility.py
  test_sec_struct.py
  test_lddt.py
  test_qsscore.py
  test_stereochemistry.py
  test_contact_score.py
  test_biounit.py
)

if (COMPOUND_LIB)
  list(APPEND OST_MOL_ALG_UNIT_TESTS test_qsscoring.py
                                     test_nonstandard.py
                                     test_chain_mapping.py
                                     test_ligand_scoring.py)
endif()

ost_unittest(MODULE mol_alg SOURCES "${OST_MOL_ALG_UNIT_TESTS}" LINK ost_io)
