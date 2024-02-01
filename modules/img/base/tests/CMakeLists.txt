set(OST_IMG_BASE_UNIT_TESTS
test_data.cc
test_extent.cc
test_image.cc
test_image_state.cc
test_image_state_visitor.cc
test_index.cc
test_point.cc
test_sampling.cc
test_domains.cc
test_size.cc
test_value_holder.cc
tests.cc
)


ost_unittest(MODULE img
             SOURCES "${OST_IMG_BASE_UNIT_TESTS}"
             LINK ost_img_alg)

