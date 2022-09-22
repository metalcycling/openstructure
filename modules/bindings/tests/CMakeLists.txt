set(OST_BINDINGS_UNIT_TESTS
  test_msms.py
  test_clustalw.py
  test_hhblits.py
  test_blast.py
  test_kclust.py
  test_naccess.py
  test_cadscore.py
  test_tmtools.py
  test_ialign.py
  test_lga.py
  test_hbplus.py
  test_dockq.py
)

ost_unittest(MODULE bindings 
             SOURCES "${OST_BINDINGS_UNIT_TESTS}")

