set(OST_SEQ_ALG_UNIT_TESTS
  test_distance_analysis.cc
  test_merge_pairwise_alignments.cc
  test_sequence_identity.cc
  tests.cc
  test_renumber.py
  test_local_align.py
  test_global_align.py
  test_semiglobal_align.py
  test_weight_matrix.py
  test_alignment_from_chain_view.py
  test_aaindex.py
)

if (COMPOUND_LIB)
  list(APPEND OST_SEQ_ALG_UNIT_TESTS test_aligntoseqres.py)
endif()

ost_unittest(MODULE seq_alg SOURCES "${OST_SEQ_ALG_UNIT_TESTS}"
             LINK ost_mol ost_io)

