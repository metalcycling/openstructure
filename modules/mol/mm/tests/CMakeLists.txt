set(OST_MOL_MM_UNIT_TESTS
  test_interaction.cc
  test_block.cc
  test_block_modifiers.cc
  test_topology.cc
  test_forcefield.cc
  test_simulation.cc
  tests.cc
  test_antechamber.py
)

ost_unittest(MODULE mol_mm SOURCES "${OST_MOL_MM_UNIT_TESTS}")

