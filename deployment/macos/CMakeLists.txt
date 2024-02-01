#-------------------------------------------------------------------------------
#  Author:  Andreas Schenk
#
#  This file contains OSX specific CPack configuration.
#
#-------------------------------------------------------------------------------


set(CPACK_GENERATOR Bundle)
set(CPACK_BUNDLE_NAME DNG)
set(CPACK_BUNDLE_PLIST ${PROJECT_SOURCE_DIR}/deployment/macos/Info.plist)
set (CPACK_PACKAGE_ICON ${PROJECT_SOURCE_DIR}/graphics/icon.icns)
set (CPACK_BUNDLE_ICON ${PROJECT_SOURCE_DIR}/graphics/icon.icns)
set(CPACK_PACKAGE_INSTALL_DIRECTORY "openstructure ${OST_VERSION_STRING}" )
install(CODE "
  if(NOT ${CMAKE_INSTALL_PREFIX} STREQUAL \${CMAKE_INSTALL_PREFIX})
    set(LIB_DIR ${LIB_DIR}) #pass libdir variable in install code block
    set(LIBEXEC_PATH ${LIBEXEC_PATH})#pass libexec_path variable in install code block
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake_support)
    include(CopyDependencies)
  
    set(OST_BINARIES \"\${CMAKE_INSTALL_PREFIX}/bin/lddt\"
                     \"\${CMAKE_INSTALL_PREFIX}/bin/chemdict_tool\"
                     \"\${CMAKE_INSTALL_PREFIX}/libexec/openstructure/gosty\")
  
    #copy python into bundle if a version different from the system version is used
    if(NOT ${Python_EXECUTABLE} MATCHES /usr/bin/)
      copy_python(${Python_INCLUDE_DIRS} ${Python_VERSION} new_python_binary)
      list(APPEND OST_BINARIES \${new_python_binary})
    endif(NOT ${Python_EXECUTABLE} MATCHES /usr/bin/)  
  
    file(GLOB_RECURSE OST_PYMODS \"\${CMAKE_INSTALL_PREFIX}/${LIB_DIR}/*.so\")
    list(APPEND  OST_BINARIES \${OST_PYMODS})
    
    #copy Qt stuff
    copy_qt(${QT_LIBRARY_DIR} ${QT_PLUGINS_DIR} qt_plugins)
    list(APPEND  OST_BINARIES \${qt_plugins})
  
    write_bundle_binary(dng ../Resources/bin/dng)
    copy_dependencies(\"\${OST_BINARIES}\") 
  endif(NOT ${CMAKE_INSTALL_PREFIX} STREQUAL \${CMAKE_INSTALL_PREFIX})
")

