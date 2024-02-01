set(OST_IMG_PYMOD_SOURCES
export_data.cc
export_data_algorithm.cc
export_extent.cc
export_image_handle.cc
export_peak.cc
export_phase.cc
export_point.cc
export_progress.cc
export_raster_image.cc
export_size.cc
export_mask.cc
export_image_list.cc
export_map.cc
wrap_img.cc
)

if (NOT ENABLE_STATIC) 
  pymod(NAME img OUTPUT_DIR ost/img 
        CPP ${OST_IMG_PYMOD_SOURCES} PY __init__.py)
endif()
