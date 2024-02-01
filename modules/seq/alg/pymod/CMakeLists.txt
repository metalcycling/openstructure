set(OST_SEQ_ALG_PYMOD_SOURCES
  wrap_seq_alg.cc
)
if (NOT ENABLE_STATIC)
  pymod(NAME seq_alg OUTPUT_DIR ost/seq/alg 
        CPP ${OST_SEQ_ALG_PYMOD_SOURCES}
        PY __init__.py mat.py renumber.py aaindex.py)
endif()

copy_if_different("${CMAKE_CURRENT_SOURCE_DIR}" "${STAGE_DIR}/share/openstructure/aaindex"
                  "aaindex1" "AAINDEX1" "ost_seq_alg_pymod")
install(FILES "aaindex1" DESTINATION "share/openstructure/aaindex/")

copy_if_different("${CMAKE_CURRENT_SOURCE_DIR}" "${STAGE_DIR}/share/openstructure/aaindex"
                  "aaindex2" "AAINDEX2" "ost_seq_alg_pymod")
install(FILES "aaindex2" DESTINATION "share/openstructure/aaindex/")

copy_if_different("${CMAKE_CURRENT_SOURCE_DIR}" "${STAGE_DIR}/share/openstructure/aaindex"
                  "aaindex3" "AAINDEX3" "ost_seq_alg_pymod")
install(FILES "aaindex3" DESTINATION "share/openstructure/aaindex/")
