set(OST_GUI_PANEL_BAR_HEADERS
button_bar.hh
drop_box.hh
panel_bar.hh
panel_manager.hh
bottom_bar.hh
button_box.hh
event_button.hh
panel_widget_container.hh
panel_bar_widget_holder.hh
side_bar.hh
splitter_panel_bar.hh
tabbed_panel_bar.hh
)

set(OST_GUI_SEQUENCE_VIEWER_HEADERS
align_properties_painter.hh
base_row.hh
background_painter.hh
conservation_painter.hh
painter.hh
secstr_row.hh
seq_secstr_painter.hh
seq_selection_painter.hh
seq_text_painter.hh
sequence_delegate.hh
sequence_model.hh
sequence_row.hh
sequence_search_bar.hh
sequence_table_view.hh
sequence_viewer.hh
tick_painter.hh
title_row.hh
alignment_view_object.hh
base_view_object.hh
sequence_view_object.hh
)

set(OST_GUI_TOOLS_HEADERS
tool.hh
tool_option_group.hh
tool_options.hh
tool_bar.hh
selection_tool.hh
measure_tool.hh
tool_option.hh
tool_manager.hh
tool_options_widget.hh
manipulator.hh
tool_options_win.hh
mouse_event.hh
map_tool.hh
)
set(OST_GUI_PYTHON_SHELL_HEADERS
main_thread_runner.hh
output_redirector.hh
completer_base.hh
dir_model.hh
gutter.hh
path_completer.hh
python_context_parser.hh
python_completer.hh
python_shell.hh
python_shell_fw.hh
python_shell_widget.hh
python_tokenizer.hh
python_syntax_highlighter.hh
text_logger.hh
python_interpreter.hh
python_interpreter_worker.hh
python_shell_text_document_layout.hh
shell_history.hh
state_machine.hh
state.hh
transition.hh
transition_guard.hh
)

set(OST_GUI_SCENE_WIN_HEADERS
context_menu.hh
current_selection_node.hh
custom_part_node.hh
entity_node.hh
entity_part_node.hh
gfx_scene_node.hh
label_node.hh
root_node.hh
render_mode_node.hh
render_modes_node.hh
scene_node.hh
scene_win.hh
scene_win_model.hh
)

set(OST_GUI_MESSAGES_HEADERS
message_widget.hh
message_box_widget.hh
log_reader.hh
)

if (ENABLE_SPNAV)
include_directories(${SPNAV_INCLUDE_DIR})
set(OST_GUI_INPUT_HEADERS
spnav_input.hh
  IN_DIR input
)
endif()

set(OST_GUI_PLOT_VIEWER_HEADERS
plot_axis_base.hh
plot_axis_horizontal.hh
plot_axis_vertical.hh
plot_data.hh
plot_data_base.hh
plot_data_entry.hh
plot_data_graphics_item.hh
plot_data_graphics_item_base.hh
plot_data_info.hh
plot_data_info_fw.hh
plot_function.hh
plot_function_info.hh
plot_function_graphics_item.hh
plot_graphics_cross_item.hh
plot_graphics_diagonal_cross_item.hh
plot_graphics_ellipse_item.hh
plot_graphics_rectangle_item.hh
plot_group_graphics_item.hh
plot_group_info.hh
plot_group.hh
plot_info_base.hh
plot_legend_panel.hh
plot_ticker.hh
plot_ticker_fw.hh
plot_types.hh
)

set(OST_GUI_DATA_VIEWER_HEADERS
argand.hh
data_viewer.hh
data_viewer_panel.hh
data_viewer_panel_base.hh
draw_overlay.hh
drawing_functions.hh
fft_panel.hh
info_panel.hh
img2qt.hh
mask_overlay.hh
null_overlay.hh
overlay.hh
overlay_base.hh
overlay_base_fw.hh
overlay_manager.hh
overlay_manager_fw.hh
overlay_manager_gui.hh
overlay_manager_observer.hh
pointlist_overlay.hh
pointlist_overlay_base.hh
strategies.hh
strategies_fw.hh
viewer_normalizer.hh
IN_DIR data_viewer
)

set(OST_GUI_HEADERS
file_browser.hh
file_type_dialog.hh
file_viewer.hh
gl_canvas.hh
gl_win.hh
scene_menu.hh
gosty_app.hh
admin.hh
module_config.hh
main.hh
main_area.hh
perspective.hh
scene_selection.hh
main_window.hh
dock_widget.hh
dock_widget.hh
thin_splitter.hh
widget.hh
widget_geom_handler.hh
widget_pool.hh
widget_registry.hh
widget_state_saver.hh
file_loader.hh
)

set(OST_GUI_SOURCES
thin_splitter.cc
perspective.cc
widget_registry.cc
file_browser.cc
file_type_dialog.cc
file_viewer.cc
main.cc
gl_canvas.cc
dock_widget.cc
main_window.cc
gl_win.cc
scene_menu.cc
widget.cc
admin.cc
widget_pool.cc
messages/message_widget.cc
messages/message_box_widget.cc
messages/log_reader.cc
sequence_viewer/align_properties_painter.cc
sequence_viewer/base_row.cc
sequence_viewer/background_painter.cc
sequence_viewer/conservation_painter.cc
sequence_viewer/secstr_row.cc
sequence_viewer/seq_secstr_painter.cc
sequence_viewer/seq_selection_painter.cc
sequence_viewer/seq_text_painter.cc
sequence_viewer/sequence_delegate.cc
sequence_viewer/sequence_model.cc
sequence_viewer/sequence_row.cc
sequence_viewer/sequence_search_bar.cc
sequence_viewer/sequence_table_view.cc
sequence_viewer/sequence_viewer.cc
sequence_viewer/tick_painter.cc
sequence_viewer/title_row.cc
sequence_viewer/alignment_view_object.cc
sequence_viewer/base_view_object.cc
sequence_viewer/sequence_view_object.cc
gosty_app.cc
main_area.cc
widget_state_saver.cc
scene_selection.cc
widget_geom_handler.cc
file_loader.cc
panels/button_bar.cc
panels/drop_box.cc
panels/panel_bar.cc
panels/bottom_bar.cc
panels/button_box.cc
panels/event_button.cc
panels/panel_manager.cc
panels/panel_widget_container.cc
panels/panel_bar_widget_holder.cc
panels/side_bar.cc
panels/splitter_panel_bar.cc
panels/tabbed_panel_bar.cc
tools/tool_option_group.cc
tools/tool_options.cc
tools/tool_option.cc
tools/tool_manager.cc
tools/tool_options_win.cc
tools/tool.cc
tools/mouse_event.cc
tools/manipulator.cc
tools/tool_options_widget.cc
tools/tool_bar.cc
tools/selection_tool.cc
tools/measure_tool.cc
python_shell/gutter.cc
python_shell/main_thread_runner.cc
python_shell/output_redirector.cc
python_shell/completer_base.cc
python_shell/dir_model.cc
python_shell/path_completer.cc
python_shell/python_context_parser.cc
python_shell/python_completer.cc
python_shell/python_interpreter.cc
python_shell/python_interpreter_worker.cc
python_shell/python_namespace_tree_model.cc
python_shell/python_namespace_tree_item.cc
python_shell/python_shell.cc
python_shell/python_shell_widget.cc
python_shell/python_shell_text_document_layout.cc
python_shell/python_tokenizer.cc
python_shell/python_syntax_highlighter.cc
python_shell/shell_history.cc
python_shell/string_literal_positions.cc
python_shell/text_logger.cc
python_shell/state_machine.cc
python_shell/state.cc
python_shell/transition.cc
python_shell/transition_guard.cc
scene_win/context_menu.cc
scene_win/current_selection_node.cc
scene_win/custom_part_node.cc
scene_win/scene_node.cc
scene_win/label_node.cc
scene_win/entity_node.cc
scene_win/entity_part_node.cc
scene_win/gfx_scene_node.cc
scene_win/render_mode_node.cc
scene_win/render_modes_node.cc
scene_win/root_node.cc
scene_win/scene_win.cc
scene_win/scene_win_model.cc
tools/map_tool.cc
data_viewer/argand.cc
data_viewer/data_viewer.cc
data_viewer/data_viewer_panel.cc
data_viewer/data_viewer_panel_base.cc
data_viewer/draw_overlay.cc
data_viewer/drawing_functions.cc
data_viewer/fft_panel.cc
data_viewer/info_panel.cc
data_viewer/img2qt.cc
data_viewer/mask_overlay.cc
data_viewer/null_overlay.cc
data_viewer/overlay_base.cc
data_viewer/overlay_manager.cc
data_viewer/overlay_manager_gui.cc
data_viewer/pointlist_overlay.cc
data_viewer/pointlist_overlay_base.cc
data_viewer/strategies.cc
data_viewer/viewer_normalizer.cc
)

set(HEADERS_TO_BE_MOCCED
file_browser.hh
file_type_dialog.hh
file_viewer.hh
gl_canvas.hh
gl_win.hh
scene_menu.hh
gosty_app.hh
main.hh
main_area.hh
perspective.hh
dock_widget.hh
main_window.hh
scene_selection.hh
widget.hh
widget_geom_handler.hh
widget_pool.hh
messages/message_widget.hh
messages/message_box_widget.hh
messages/log_reader.hh
sequence_viewer/align_properties_painter.hh
sequence_viewer/background_painter.hh
sequence_viewer/base_row.hh
sequence_viewer/conservation_painter.hh
sequence_viewer/painter.hh
sequence_viewer/secstr_row.hh
sequence_viewer/seq_secstr_painter.hh
sequence_viewer/seq_selection_painter.hh
sequence_viewer/seq_text_painter.hh
sequence_viewer/sequence_delegate.hh
sequence_viewer/sequence_model.hh
sequence_viewer/sequence_row.hh
sequence_viewer/sequence_search_bar.hh
sequence_viewer/sequence_table_view.hh
sequence_viewer/sequence_viewer.hh
sequence_viewer/tick_painter.hh
sequence_viewer/title_row.hh
sequence_viewer/alignment_view_object.hh
sequence_viewer/base_view_object.hh
sequence_viewer/sequence_view_object.hh
python_shell/main_thread_runner.hh
python_shell/completer_base.hh
python_shell/gutter.hh
python_shell/output_redirector.hh
python_shell/path_completer.hh
python_shell/python_completer.hh
python_shell/python_interpreter.hh
python_shell/python_interpreter_worker.hh
python_shell/python_namespace_tree_model.hh
python_shell/python_shell.hh
python_shell/python_shell_widget.hh
python_shell/text_logger.hh
python_shell/state_machine.hh
python_shell/state.hh
python_shell/transition.hh
python_shell/transition_guard.hh
panels/button_bar.hh
panels/drop_box.hh
panels/panel_bar.hh
panels/panel_manager.hh
panels/bottom_bar.hh
panels/button_box.hh
panels/event_button.hh
panels/panel_widget_container.hh
panels/panel_bar_widget_holder.hh
panels/side_bar.hh
panels/splitter_panel_bar.hh
panels/tabbed_panel_bar.hh
scene_win/context_menu.hh
scene_win/current_selection_node.hh
scene_win/custom_part_node.hh
scene_win/gfx_scene_node.hh
scene_win/label_node.hh
scene_win/entity_node.hh
scene_win/entity_part_node.hh
scene_win/render_mode_node.hh
scene_win/render_modes_node.hh
scene_win/root_node.hh
scene_win/scene_node.hh
scene_win/scene_win.hh
scene_win/scene_win_model.hh
tools/tool.hh
tools/tool_option_group.hh
tools/tool_options.hh
tools/tool_bar.hh
tools/tool_option.hh
tools/tool_manager.hh
tools/tool_options_widget.hh
tools/measure_tool.hh
tools/tool_options_win.hh
data_viewer/argand.hh
data_viewer/fft_panel.hh
data_viewer/info_panel.hh
data_viewer/data_viewer.hh
data_viewer/overlay_manager_gui.hh
data_viewer/data_viewer_panel_base.hh
data_viewer/null_overlay.hh
data_viewer/overlay_base.hh
data_viewer/pointlist_overlay_base.hh
)

if (ENABLE_SPNAV)
  set(spnav_enabled 1)
  list(APPEND OST_GUI_SOURCES
      input/spnav_input.cc
  )
  list(APPEND HEADERS_TO_BE_MOCCED
      input/spnav_input.hh
  )
else()
  set(spnav_enabled 0)
endif()

set(config_hh_generator "CMake")
set(MODULE_CONFIG_HH_FILE "${CMAKE_CURRENT_SOURCE_DIR}/module_config.hh")
configure_file(module_config.hh.in ${MODULE_CONFIG_HH_FILE})
# the resulting file, module_config.hh, is already included in OST_GUI_HEADERS

if (APPLE)
  set(ADDITIONAL_LIBRARIES "-framework Security")
endif()
set(QT_USE_QTOPENGL 1)
set(QT_USE_QTNETWORK 1)

qt5_wrap_cpp(OST_GUI_MOCS "${HEADERS_TO_BE_MOCCED}" 
             OPTIONS -DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED
                     -DBOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)
module(NAME gui SOURCES ${OST_GUI_MOCS} ${OST_GUI_SOURCES}
       HEADERS ${OST_GUI_TOOLS_HEADERS} IN_DIR tools
               ${OST_GUI_SEQUENCE_VIEWER_HEADERS} IN_DIR sequence_viewer
               ${OST_GUI_PYTHON_SHELL_HEADERS} IN_DIR python_shell
               ${OST_GUI_PANEL_BAR_HEADERS} IN_DIR panels
               ${OST_GUI_SCENE_WIN_HEADERS} IN_DIR scene_win
               ${OST_GUI_MESSAGES_HEADERS} IN_DIR messages
               ${OST_GUI_INPUT_HEADERS}
               ${OST_GUI_DATA_VIEWER_HEADERS}
               ${OST_GUI_HEADERS}
       DEPENDS_ON ost_gfx ost_io ost_mol_alg ost_seq_alg
       LINK ${QT_LIBRARIES}
            ${Python_LIBRARIES} 
            ${BOOST_PYTHON_LIBRARIES}
            ${SPNAV_LIBRARIES}
        NO_STATIC)


   
if (ADDITIONAL_LIBRARIES)
  target_link_libraries(ost_gui "${ADDITIONAL_LIBRARIES}")
endif()

include_directories(${Python_INCLUDE_DIRS})
qt5_add_resources(OST_QT_RESOURCE dngr.qrc)
qt5_wrap_cpp(OST_GOSTY_MOC "gosty.hh"
             OPTIONS -DBOOST_TT_HAS_OPERATOR_HPP_INCLUDED
                     -DBOOST_NO_TEMPLATE_PARTIAL_SPECIALIZATION)

set(LINK LINK ${QT_LIBRARIES} ${BOOST_PROGRAM_OPTIONS_LIBRARIES})

executable_libexec(NAME gosty SOURCES gosty.cc ${OST_GOSTY_MOC}
                   ${OST_QT_RESOURCE} DEPENDS_ON ost_gui ${LINK})
           
file(GLOB MOC_CXX_FILES moc_*.cxx)

if (MOC_CXX_FILES)
  set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
                           "${MOC_CXX_FILES}")
endif()

