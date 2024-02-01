set(IPLT_ALG_UNIT_TESTS
test_clear.cc
test_conj.cc
test_convolute.cc
test_discrete_shrink.cc
test_fft.cc
test_fill.cc
test_filter.cc

test_histogram.cc
test_mirror.cc
test_negate.cc
test_power_spectrum.cc
test_shift.cc
test_stat.cc
test_transform.cc
test_normalizer.cc
tests.cc
)

ost_unittest(MODULE img_alg
             SOURCES "${IPLT_ALG_UNIT_TESTS}")
