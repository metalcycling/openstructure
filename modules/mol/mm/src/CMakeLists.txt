set(OST_MOL_MM_HEADERS
  forcefield.hh
  ff_reader.hh
  buildingblock.hh
  block_modifiers.hh
  gromacs_block_modifiers.hh
  heuristic_block_modifiers.hh
  interaction.hh
  settings.hh
  simulation.hh 
  modeller.hh
  observer.hh
  state_extractor.hh
  system_creator.hh
  topology_creator.hh
  index.hh
  topology.hh
  steep.hh
)

set(OST_MOL_MM_SOURCES
  forcefield.cc
  ff_reader.cc
  buildingblock.cc
  gromacs_block_modifiers.cc
  heuristic_block_modifiers.cc
  interaction.cc
  simulation.cc
  modeller.cc
  state_extractor.cc
  observer.cc
  system_creator.cc
  topology_creator.cc
  topology.cc
  steep.cc
)

# create settings.hh as configurational header, needed to set the plugins path
set(config_hh_generator "CMake")
set(SETTINGS_HH_FILE "${CMAKE_CURRENT_SOURCE_DIR}/settings.hh")
configure_file(settings.hh.in ${SETTINGS_HH_FILE})


set(MOL_MM_DEPS ost_mol ost_io)


module(NAME mol_mm SOURCES ${OST_MOL_MM_SOURCES}
       HEADERS ${OST_MOL_MM_HEADERS}
       HEADER_OUTPUT_DIR ost/mol/mm
       DEPENDS_ON ${MOL_MM_DEPS}
       LINK ${BOOST_PROGRAM_OPTIONS})

target_link_libraries(ost_mol_mm ${OPEN_MM_LIBRARIES})

copy_if_different("${CMAKE_CURRENT_SOURCE_DIR}" "${STAGE_DIR}/share/openstructure/forcefields"
                  "AMBER03.dat" "AMBER_03_FORCEFIELD"
                  "ost_mol_mm")
install(FILES "AMBER03.dat" DESTINATION "share/openstructure/forcefields/")

copy_if_different("${CMAKE_CURRENT_SOURCE_DIR}" "${STAGE_DIR}/share/openstructure/forcefields"
                  "CHARMM27.dat" "CHARMM_27_FORCEFIELD"
                  "ost_mol_mm")
install(FILES "CHARMM27.dat" DESTINATION "share/openstructure/forcefields/")

