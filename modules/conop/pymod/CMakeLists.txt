set(OST_CONOP_PYMOD_SOURCES
  wrap_conop.cc
  #  export_builder.cc
  export_compound.cc
  export_processor.cc
  export_heuristic.cc
  export_amino_acids.cc
  export_conop.cc
  export_diag.cc
  export_rule_based.cc
  export_ring_finder.cc
)

if (NOT ENABLE_STATIC)
  pymod(NAME conop CPP ${OST_CONOP_PYMOD_SOURCES} PY __init__.py cleanup.py)
endif()
