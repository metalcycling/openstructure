set(OST_DB_UNIT_TESTS
  test_linear_db.py
)

ost_unittest(MODULE db 
             SOURCES "${OST_DB_UNIT_TESTS}")

