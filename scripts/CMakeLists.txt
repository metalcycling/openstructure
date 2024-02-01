set(SUBST_DICT BUILD_TYPE=${CMAKE_BUILD_TYPE} Python_EXECUTABLE=${Python_EXECUTABLE}
Python_VERSION=${Python_VERSION} LIBDIR=${LIB_DIR} LIBEXEC_PATH=${LIBEXEC_PATH}
PYTHON_MODULE_PATH=${PYTHON_MODULE_PATH})
message("${BUILD_TYPE}")
add_custom_target(ost_scripts ALL)

# may god forgive us
if (ENABLE_GUI)
  add_dependencies(ost_scripts _ost_gui)
elseif(ENABLE_GFX)
  add_dependencies(ost_scripts _ost_gfx)
else()
  add_dependencies(ost_scripts _ost_io)
endif()

if (NOT ENABLE_STATIC)
  if(WIN32)
    script(NAME  dng.bat INPUT dng.bat.in SUBSTITUTE ${SUBST_DICT}
           TARGET ost_scripts)
    script(NAME ost.bat INPUT ost.bat.in SUBSTITUTE ${SUBST_DICT}
           TARGET ost_scripts)
  else()
    script(NAME ost_config INPUT ost_config.in SUBSTITUTE ${SUBST_DICT}
           TARGET ost_scripts OUTPUT_DIR ${LIBEXEC_PATH})
    script(NAME ost INPUT ost.in SUBSTITUTE ${SUBST_DICT} 
           TARGET ost_scripts)
    if (ENABLE_GUI)
      script(NAME dng INPUT dng.in SUBSTITUTE ${SUBST_DICT}
             TARGET ost_scripts)
    endif()
  endif()

  script(NAME ost_startup.py INPUT ost_startup.py.in SUBSTITUTE ${SUBST_DICT} 
         TARGET ost_scripts OUTPUT_DIR "${LIB_DIR}/${PYTHON_MODULE_PATH}/ost" )    

  set(INIT_SCRIPTS 
  )
  set(SPD "${LIB_DIR}/${PYTHON_MODULE_PATH}/ost")
  if (ENABLE_GUI)
    list(APPEND INIT_SCRIPTS dng_startup.py)
  endif()
  copy_if_different("${CMAKE_CURRENT_SOURCE_DIR}" "${STAGE_DIR}/${SPD}" 
                  "${INIT_SCRIPTS}" "python init scripts" ost_scripts)

  install(FILES ${INIT_SCRIPTS} DESTINATION "${SPD}")

endif()
