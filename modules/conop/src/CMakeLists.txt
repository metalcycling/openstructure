set(OST_CONOP_HEADERS
conop.hh
processor.hh
amino_acids.hh
diag.hh
model_check.hh
heuristic.hh
compound.hh
compound_lib.hh
module_config.hh
rule_based.hh
minimal_compound_lib.hh
compound_lib_base.hh
ring_finder.hh
)

set(OST_CONOP_SOURCES
processor.cc
amino_acids.cc
conop.cc
minimal_compound_lib.cc
compound_lib_base.cc
heuristic.cc
diag.cc
rule_based.cc
model_check.cc
compound.cc
compound_lib.cc
ring_finder.cc
)

module(NAME conop SOURCES ${OST_CONOP_SOURCES}
       HEADERS ${OST_CONOP_HEADERS} DEPENDS_ON ost_mol ost_geom ost_db)


if (WIN32)
  executable(NAME chemdict_tool SOURCES chemdict_tool.cc ../../io/src/zlib.cpp ../../io/src/gzip.cpp DEPENDS_ON ost_io STATIC)
else(WIN32)
  executable(NAME chemdict_tool SOURCES chemdict_tool.cc DEPENDS_ON ost_io STATIC)
endif(WIN32)

if (COMPOUND_LIB)
  if (EXISTS "${COMPOUND_LIB}")
    # this ensures that file is called 'compounds.chemlib' (as required)
    set(TO "${STAGE_DIR}/share/openstructure/compounds.chemlib")
    file(MAKE_DIRECTORY "${STAGE_DIR}/share/openstructure")
    add_custom_command(TARGET "ost_conop" PRE_BUILD
          DEPENDS "${COMPOUND_LIB}"
          COMMAND ${CMAKE_COMMAND} -E copy_if_different ${COMPOUND_LIB} ${TO})
    install(FILES ${TO} DESTINATION "share/openstructure")
  else()
    message(FATAL_ERROR "${COMPOUND_LIB} does not exist. \n"
            "Pass -DCOMPOUND_LIB=OFF if you would like to build "
            "OpenStructure without a compound library")
  endif()

endif()
