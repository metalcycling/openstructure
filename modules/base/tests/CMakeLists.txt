set(OST_BASE_UNIT_TESTS
  test_generic_property.cc
  test_string_ref.cc
  test_pod_vector.cc
  test_stutil.py
  test_table.py
  test_log.py
  tests.cc
)

ost_unittest(MODULE base 
             SOURCES "${OST_BASE_UNIT_TESTS}"
             LINK ost_mol)
