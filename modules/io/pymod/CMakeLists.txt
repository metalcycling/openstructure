set(OST_IO_PYMOD_SOURCES
  wrap_io.cc
  export_pdb_io.cc
  export_mmcif_io.cc
  export_omf_io.cc
  export_map_io.cc
)

set(OST_IO_PYMOD_MODULES
  __init__.py
  remote.py
  repository.py
)

if (NOT ENABLE_STATIC)
  pymod(NAME io CPP ${OST_IO_PYMOD_SOURCES} PY ${OST_IO_PYMOD_MODULES})
  if (STATIC_LIBRARIES)
    message("adding zlib libraries")
    target_link_libraries(_ost_io ${ZLIB_LIBRARIES})
  endif()
endif()

