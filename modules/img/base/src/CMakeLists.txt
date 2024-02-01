add_subdirectory(raster_image)
add_subdirectory(image_state)


set(OST_IMG_SOURCES 
circle_mask.cc
composite_mask.cc
data_algorithm.cc
data_observer.cc
data.cc
extent_iterator.cc
extent_mask.cc
extent.cc
image_factory.cc
image_handle.cc
image_impl.cc
image_list.cc
inverted_mask.cc
map.cc
mask_base.cc
mask_op.cc
mask.cc
observable.cc
paste_impl.cc
peak.cc
phase.cc
physical_units.cc
pixel_sampling.cc
point_list.cc
point.cc
polygon_mask.cc
progress.cc
size.cc
spherical_mask.cc
)

set(OST_IMG_HEADERS
algorithm.hh
circle_mask.hh
composite_mask.hh
data_algorithm.hh
data_observer.hh
data_types.hh
data.hh
debug.hh
dllexport.hh
extent_iterator.hh
extent_mask.hh
extent.hh
image_factory.hh
image_fw.hh
image_handle_fw.hh
image_handle.hh
image_impl.hh
image_list.hh
image_state_fw.hh
image_state.hh
image.hh
inverted_mask.hh
map.hh
mask_base_fw.hh
mask_base.hh
mask_op.hh
mask_visitor.hh
mask.hh
module_config.hh
normalizer_fw.hh
normalizer_impl.hh
null_data.hh
observable.hh
paste_impl.hh
peak.hh
phase.hh
physical_units.hh
pixel_sampling.hh
point_list.hh
point.hh
polygon_mask.hh
progress.hh
size.hh
spherical_mask.hh
util.cc
util.hh
value_util.hh
vecmat.hh
)

if (ENABLE_INFO)
  list(APPEND OST_IMG_SOURCES mask_info_convert.cc)
  list(APPEND OST_IMG_HEADERS mask_info_convert.hh)
  set(INFO_DEPS ost_info)
endif()

foreach(fname ${OST_IMG_IMAGE_STATE_SOURCES})
  set(OST_IMG_SOURCES ${OST_IMG_SOURCES} image_state/${fname})
endforeach(fname ${OST_IMG_IMAGE_STATE_SOURCES})

foreach(fname ${OST_IMG_RASTER_IMAGE_SOURCES})
  set(OST_IMG_SOURCES ${OST_IMG_SOURCES} raster_image/${fname})
endforeach(fname ${OST_IMG_RASTER_IMAGE_SOURCES})

module(NAME img SOURCES "${OST_IMG_SOURCES}"
       HEADERS ${OST_IMG_RASTER_IMAGE_HEADERS} IN_DIR raster_image
                ${OST_IMG_IMAGE_STATE_HEADERS} IN_DIR image_state
                ${OST_IMG_HEADERS}
       DEPENDS_ON ost_geom ost_base ${INFO_DEPS})
