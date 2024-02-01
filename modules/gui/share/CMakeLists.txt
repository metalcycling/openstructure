
set(GUI_ICONS
  icons/add_icon.png
  icons/close_icon.png
  icons/critical_icon.png
  icons/delete_icon.png
  icons/distance_icon.png
  icons/find_icon.png
  icons/information_icon.png
  icons/map_icon.png
  icons/menubar_icon.png
  icons/rigid_body_manipulator_icon.png
  icons/selection_icon.png
  icons/show_sidebar_icon.png
  icons/site_icon.png
  icons/split_icon.png
  icons/warning_icon.png
  icons/refresh.png
)

set(GUI_IMAGES
  images/logo-small.png
)

add_custom_target(scene_icons ALL)
copy_if_different("${CMAKE_CURRENT_SOURCE_DIR}" "${STAGE_DIR}/share/openstructure/gui/icons" 
                  "${GUI_ICONS}" "ICONS_FOR_GUI" scene_icons)
install(FILES ${GUI_ICONS} DESTINATION "share/openstructure/gui/icons")

add_custom_target(gui_images ALL)
copy_if_different("${CMAKE_CURRENT_SOURCE_DIR}" "${STAGE_DIR}/share/openstructure/gui/images" 
                  "${GUI_IMAGES}" "IMAGES_FOR_GUI" gui_images)
install(FILES ${GUI_IMAGES} DESTINATION "share/openstructure/gui/images")