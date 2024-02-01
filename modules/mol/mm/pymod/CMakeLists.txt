set(OST_MOL_MM_PYMOD_SOURCES
  export_forcefield.cc
  export_settings.cc
  export_simulation.cc
  export_openmm.cc
  export_observers.cc
  export_ff_reader.cc
  export_interaction.cc
  export_buildingblock.cc
  export_block_modifiers.cc
  export_topology.cc
  export_modeller.cc
  wrap_mol_mm.cc
)

set(OST_MOL_MM_PYMOD_MODULES
  __init__.py
  antechamber.py
)

if (NOT ENABLE_STATIC)
  pymod(NAME mol_mm OUTPUT_DIR ost/mol/mm CPP ${OST_MOL_MM_PYMOD_SOURCES}
        PY ${OST_MOL_MM_PYMOD_MODULES})
endif()
