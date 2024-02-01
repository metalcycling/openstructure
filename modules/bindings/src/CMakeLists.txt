if (COMPILE_TMTOOLS)
  enable_language(Fortran)
  executable(NAME tmalign SOURCES tmalign.f)
  executable(NAME tmscore SOURCES tmscore.f)
endif()

set(OST_BINDINGS_SOURCES 
wrap_tmalign.cc)

set(OST_BINDINGS_HEADERS
wrap_tmalign.hh)

module(NAME bindings SOURCES ${OST_BINDINGS_SOURCES}
       HEADERS ${OST_BINDINGS_HEADERS} HEADER_OUTPUT_DIR ost/bindings
       DEPENDS_ON ost_geom ost_mol ost_seq)
