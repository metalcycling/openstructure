set(OST_GFX_PYMOD_SOURCES
  wrap_gfx.cc
  export_scene.cc
  export_gfx_node.cc
  export_symmetry_node.cc
  export_gfx_obj.cc
  export_entity.cc
  export_surface.cc
  export_primitives.cc
  export_primlist.cc
  export_scene_observer.cc
  export_render_options.cc
  export_color.cc
  export_gradient.cc
  export_color_ops.cc
  export_glwin_base.cc
  export_exporter.cc
  export_bitmap.cc
  export_map.cc
)

pymod(NAME gfx CPP ${OST_GFX_PYMOD_SOURCES} PY __init__.py py_gfx_obj.py)

set(GRADIENT_FILE
  gradients.xml
)
copy_if_different("${CMAKE_CURRENT_SOURCE_DIR}" "${STAGE_DIR}/share/openstructure/scene" 
                  "${GRADIENT_FILE}" "PRESET GRADIENTS" _ost_gfx)
                  
install(FILES ${GRADIENT_FILE} DESTINATION "share/openstructure/scene")
