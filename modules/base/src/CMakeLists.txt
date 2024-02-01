set(OST_BASE_SOURCES 
generic_property.cc
log.cc
log_sink.cc
profile.cc
units.cc
string_ref.cc
platform.cc
message.cc
test_utils/compare_files.cc
boost_filesystem_helper.cc
)

set(OST_BASE_HEADERS
generic_property.hh
integrity_error.hh
invalid_handle.hh
ptr_observer.hh
log.hh
log_sink.hh
profile.hh
dyn_cast.hh
stdint.hh
stdint_msc.hh
range.hh
module_config.hh
units.hh
message.hh
platform.hh
string_ref.hh
pod_vector.hh
fixed_string.hh
tri_matrix.hh
boost_filesystem_helper.hh
paged_array.hh
)

set(OST_EXPORT_HELPERS
generic_property_def.hh 
pair_to_tuple_conv.hh 
vec_to_list_conv.hh
)
module(NAME base SOURCES ${OST_BASE_SOURCES} 
       HEADERS ${OST_EXPORT_HELPERS} IN_DIR export_helper 
       compare_files.hh IN_DIR test_utils ${OST_BASE_HEADERS}
       DEPENDS_ON ost_geom
       HEADER_OUTPUT_DIR ost)

