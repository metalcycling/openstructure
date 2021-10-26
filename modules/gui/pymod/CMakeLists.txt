set(OST_GUI_PYMOD_SOURCES
  wrap_gui.cc
  export_alignment_view.cc
  export_message_widget.cc  
  export_gl_win.cc
  export_tool.cc
  export_py_shell.cc
  export_gosty.cc
  export_scene_win.cc
  export_sequence_viewer.cc
  export_perspective.cc
  export_sip_handler.cc
  export_scene_selection.cc
  export_main_area.cc
  export_panels.cc
  export_menu_bar.cc
  export_file_loader.cc
  export_file_viewer.cc
  export_widget.cc
  export_gl_canvas.cc
  data_viewer_proxy.cc
  export_data_viewer.cc
  export_overlay.cc
  export_overlay_manager.cc
)
set(OST_GUI_SCENE_PYMOD_MODULES
__init__.py
color_options_widget.py
color_select_widget.py
combo_options_widget.py
cpk_widget.py
custom_widget.py
gradient_editor_widget.py
gradient_info_handler.py
gradient_list_model.py
gradient_preset_widget.py
hsc_widget.py
immutable_gradient_info_handler.py
immutable_preset_info_handler.py
init_inspector.py
inspector_widget.py
map_level_widget.py
preset.py
remote.py
preset_editor_list_model.py
preset_editor_widget.py
preset_info_handler.py
preset_list_model.py
preset_widget.py
render_mode_widget.py
render_op.py
render_options_widget.py
scene_observer_impl.py
scene_selection_helper.py
simple_widget.py
sline_widget.py
toolbar_options_widget.py
trace_widget.py
tube_widget.py
uniform_color_widget.py
visibility_op.py
line_trace_widget.py
wireframe_widget.py
query_editor.py
)

if (ENABLE_SPNAV)
  list( APPEND OST_GUI_PYMOD_SOURCES
  export_input.cc
  )
endif()

set(OST_GUI_PYMOD_MODULES
  __init__.py
  init_context_menu.py
  init_menubar.py
  init_spacenav.py
  init_splash.py
  trajectory_viewer.py
  helpwidget.py
  table.py
  gosty_startup.py
)

set(OST_GUI_PYMOD_DNG_MODULES
  __init__.py
  termuse.py
  superpositiondialog.py
  init.py
  menu.py
)

pymod(NAME gui CPP ${OST_GUI_PYMOD_SOURCES} 
      PY ${OST_GUI_SCENE_PYMOD_MODULES} IN_DIR scene
         ${OST_GUI_PYMOD_DNG_MODULES} IN_DIR dng
         ${OST_GUI_PYMOD_MODULES})

set(PRESET_FILES
  scene/presets.xml
)
copy_if_different("${CMAKE_CURRENT_SOURCE_DIR}" "${STAGE_DIR}/share/openstructure/scene" 
                  "${PRESET_FILES}" "IMMUTABLE_FILES" _ost_gui)
install(FILES ${PRESET_FILES} DESTINATION "share/openstructure/scene")

set(INSPECTOR_ICONS
  scene/icons/color_icon.png
  scene/icons/preset_icon.png  
  scene/icons/render_icon.png
  scene/icons/tool_icon.png
)
copy_if_different("${CMAKE_CURRENT_SOURCE_DIR}" "${STAGE_DIR}/share/openstructure/scene/icons" 
                  "${INSPECTOR_ICONS}" "ICONS_FOR_INSPECTOR" _ost_gui)
install(FILES ${INSPECTOR_ICONS} DESTINATION "share/openstructure/scene/icons")

set(QT_USE_QTOPENGL 1)
set(QT_USE_QTNETWORK 1)

target_link_libraries(_ost_gui ${QT_LIBRARIES})
