set(OST_SEQ_IMPL_HEADERS
sequence_impl.hh
sequence_impl_fw.hh
sequence_list_impl_fw.hh
sequence_list_impl.hh
sequence_list_iterator.hh
)
set(OST_SEQ_HEADERS
sequence_list.hh
sequence_handle.hh
alignment_handle.hh
module_config.hh
sequence_op.hh
aligned_region.hh
aligned_column.hh
aligned_column_iterator.hh
invalid_sequence.hh
views_from_sequences.hh
profile_handle.hh
)

set(OST_SEQ_SOURCES
impl/sequence_impl.cc
impl/sequence_list_impl.cc
sequence_handle.cc
aligned_region.cc
aligned_column.cc
sequence_list.cc
alignment_handle.cc
sequence_op.cc
views_from_sequences.cc
profile_handle.cc
)

if (ENABLE_INFO)
  set(INFO_DEPS ost_info)
endif()


module(NAME seq SOURCES ${OST_SEQ_SOURCES} 
         HEADERS ${OST_SEQ_IMPL_HEADERS} IN_DIR impl
                 ${OST_SEQ_HEADERS}
         LINK ${BOOST_IOSTREAM_LIBRARIES}
         DEPENDS_ON ost_mol ${INFO_DEPS})