set(OST_IMG_ALG_PYMOD_SOURCES
wrap_alg.cc
export_filter.cc
export_normalizer.cc
export_transcendentals.cc
export_polar.cc
)

if (NOT ENABLE_STATIC)
  pymod(NAME img_alg OUTPUT_DIR ost/img/alg CPP ${OST_IMG_ALG_PYMOD_SOURCES}
        PY __init__.py)
endif()
