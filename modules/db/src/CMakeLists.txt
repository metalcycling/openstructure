set(OST_DB_HEADERS
sqlite_wrap.hh
sqlite_conv.hh
module_config.hh
linear_indexer.hh
binary_container.hh
extract_data_helper.hh
)

set(OST_DB_SOURCES
sqlite_wrap.cc
linear_indexer.cc
binary_container.cc
extract_data_helper.cc
)

module(NAME db SOURCES ${OST_DB_SOURCES} HEADERS ${OST_DB_HEADERS} 
       DEPENDS_ON ost_base ost_geom ost_seq LINK ${SQLITE3_LIBRARIES})
if(WIN32)
  set_target_properties(ost_db PROPERTIES LINK_FLAGS "/DEF:${CMAKE_CURRENT_SOURCE_DIR}/sqlite3.def")
  add_definitions(/DSQLITE_ENABLE_COLUMN_METADATA)
else(WIN32)
  add_definitions(-DSQLITE_OMIT_LOAD_EXTENSION)
endif(WIN32)
if (ENABLE_STATIC AND UNIX AND NOT APPLE)
  target_link_libraries(ost_db ${CMAKE_DL_LIBS})
endif()
