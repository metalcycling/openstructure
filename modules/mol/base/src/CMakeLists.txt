add_subdirectory(impl)

set(OST_MOL_SOURCES 
atom_base.cc
atom_handle.cc
atom_view.cc
bond_handle.cc
chain_base.cc
chain_handle.cc
chain_view.cc
chain_type.cc
coord_frame.cc
coord_group.cc
editor_base.cc
entity_base.cc
entity_handle.cc
entity_observer.cc
entity_property_mapper.cc
entity_view.cc
entity_visitor.cc
ics_editor.cc
not_connected_error.cc
property_id.cc
query.cc
query_error.cc
query_state.cc
residue_base.cc
residue_handle.cc
residue_view.cc
surface_builder.cc
surface_handle.cc
transfer_connectivity.cc
torsion_handle.cc
query_view_wrapper.cc
view_op.cc
coord_source.cc
in_mem_coord_source.cc
bounding_box.cc
xcs_editor.cc)

set(OST_MOL_HEADERS
atom_base.hh
atom_handle.hh
atom_view.hh
bond_handle.hh
bond_table.hh
chain_base.hh
chain_handle.hh
chain_view.hh
chain_type.hh
chem_class.hh
chem_type.hh
coord_group.hh
coord_source.hh
in_mem_coord_source.hh
coord_frame.hh
editor_base.hh
editor_type_fw.hh
mol.hh
entity_base.hh
entity_base.hh
entity_handle.hh
entity_observer.hh
entity_observer_fw.hh
entity_property_mapper.hh
entity_view.hh
entity_visitor.hh
entity_visitor_fw.hh
handle_type_fw.hh
ics_editor.hh
not_connected_error.hh
property_id.hh
query.hh
query_error.hh
query_state.hh
residue_base.hh
residue_handle.hh
residue_prop.hh
residue_view.hh
sec_structure.hh
spatial_organizer.hh
surface.hh
surface_builder.hh
builder.hh
surface_handle.hh
surface_prop.hh
torsion_handle.hh
transfer_connectivity.hh
query_view_wrapper.hh
view_op.hh
view_type_fw.hh
view_op.hh
xcs_editor.hh
bounding_box.hh
transform.hh
module_config.hh)

foreach(_impl_src ${OST_MOL_IMPL_SOURCES})
  list(APPEND OST_MOL_SOURCES impl/${_impl_src})
endforeach()

if(ENABLE_INFO)
  set (INFO_DEPS ost_info)
else()
  list(APPEND LINK ${Boost_REGEX_LIBRARY})
endif()

module(NAME mol SOURCES ${OST_MOL_SOURCES}
       HEADERS ${OST_MOL_IMPL_HEADERS} IN_DIR impl
       ${OST_MOL_HEADERS} HEADER_OUTPUT_DIR ost/mol
       DEPENDS_ON ost_geom ost_base ${INFO_DEPS} LINK ${LINK})
