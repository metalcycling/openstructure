set(OST_CONOP_UNIT_TESTS
  test_amino_acids.cc
  test_heuristic_conop.cc
  tests.cc
  test_rule_based_conop.cc 
  helper.cc
  test_processor.py
)

if (COMPOUND_LIB)
  list(APPEND OST_CONOP_UNIT_TESTS test_compound.py
                                   test_cleanup.py
                                   test_complib.py)
endif()

ost_unittest(MODULE conop
             SOURCES "${OST_CONOP_UNIT_TESTS}")
