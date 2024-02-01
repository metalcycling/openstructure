set(OST_GFX_UNIT_TESTS
  tests.cc
  test_gfx_node.cc
  test_color.cc
  test_gfx.py
  test_map_octree.cc
)

ost_unittest(MODULE gfx
             SOURCES "${OST_GFX_UNIT_TESTS}"
             LINK ost_io)

