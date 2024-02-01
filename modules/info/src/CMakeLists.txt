set(OST_INFO_SOURCES 
info_impl.cc
info_handle.cc
info_group.cc
info_item.cc
info_path.cc
info_set.cc
info_visitor.cc
geom_info_conversion.cc
info_verbose_visitor.cc)

set(OST_INFO_HEADERS
info.hh
info_error.hh
info_impl.hh
info_impl_fw.hh
info_handle.hh
info_group.hh
info_item.hh
info_fw.hh
module_config.hh
info_path.hh
info_set.hh
info_visitor.hh
info_verbose_visitor.hh
item_type.hh
geom_info_conversion.hh
item_type_cast.hh)

set(QT_USE_QTXML 1)

module(NAME info SOURCES ${OST_INFO_SOURCES}
       HEADERS ${OST_INFO_HEADERS}
       DEPENDS_ON ost_base LINK ${BOOST_REGEX_LIBRARIES} ${QT_LIBRARIES})

