#-------------------------------------------------------------------------------
#  Author:  Andreas Schenk
#
#  This file contains the Linux specific CPack configuration.
#
#-------------------------------------------------------------------------------

if(NOT CPACK_GENERATOR)
  include(CheckLSBTypes)
  if(( ${LSB_DISTRIBUTOR_ID} MATCHES "centos" AND ${LSB_RELEASE} GREATER 5.9)
   OR (${LSB_DISTRIBUTOR_ID} MATCHES "redhat" AND ${LSB_RELEASE} GREATER 5.9)
   OR (${LSB_DISTRIBUTOR_ID} MATCHES "fedora"))
    set(CPACK_GENERATOR "RPM")
  elseif(${LSB_DISTRIBUTOR_ID} MATCHES "debian")
    set(CPACK_GENERATOR "DEB")
  else()
    set(CPACK_GENERATOR "TGZ")
  endif()
endif(NOT CPACK_GENERATOR)

if(CPACK_GENERATOR MATCHES "DEB")
  set(CPACK_DEBIAN_PACKAGE_MAINTAINER "The OpenStructure developers")
elseif(CPACK_GENERATOR MATCHES "RPM")
elseif(CPACK_GENERATOR MATCHES "TGZ")
  install(CODE "
    if(NOT ${CMAKE_INSTALL_PREFIX} MATCHES \${CMAKE_INSTALL_PREFIX})
      set(LIB_DIR ${LIB_DIR})#pass libdir variable in install code block
      set(LIBEXEC_PATH ${LIBEXEC_PATH})#pass libexec_path variable in install code block
      set(_UBUNTU_LAYOUT ${_UBUNTU_LAYOUT})#pass ubuntu_layout variable in install code block
      set(ARCH ${ARCH})#pass arch variable in install code block
      if(NOT ${_USE_RPATH} STREQUAL ON)
        message(\"======================================================================\")
        message(\"|          Warning, no rpaths will be embedded upon install.         |\")
        message(\"|        Dependency resolution during packaging will not work.       |\")
        message(\"| Add -DUSE-RPATH=1 when running cmake to switch on rpath embedding. |\")
        message(\"======================================================================\")
      endif(NOT ${_USE_RPATH} STREQUAL ON)
      set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/ost/cmake_support)
      include(CopyDependencies)
      set(OST_BINARIES \"\${CMAKE_INSTALL_PREFIX}/bin/lddt\"
                       \"\${CMAKE_INSTALL_PREFIX}/bin/chemdict_tool\"
                       \"\${CMAKE_INSTALL_PREFIX}/bin/molck\"
                       \"\${CMAKE_INSTALL_PREFIX}/${LIBEXEC_PATH}/gosty\")

      copy_python(${Python_INCLUDE_DIRS} ${Python_VERSION} new_python_binary)
      list(APPEND OST_BINARIES \${new_python_binary})

      file(GLOB_RECURSE OST_PYMODS \"\${CMAKE_INSTALL_PREFIX}/${LIB_DIR}/*.so\")
      list(APPEND  OST_BINARIES \${OST_PYMODS})

      #copy Qt stuff
      copy_qt(${QT_LIBRARY_DIR} ${QT_PLUGINS_DIR} qt_plugins)
      list(APPEND  OST_BINARIES \${qt_plugins})

      copy_dependencies(\"\${OST_BINARIES}\")

      #copy examples into bundle
      file (COPY ${CMAKE_SOURCE_DIR}/examples DESTINATION \"\${CMAKE_INSTALL_PREFIX}/share/openstructure/\") 
    endif(NOT ${CMAKE_INSTALL_PREFIX} MATCHES \${CMAKE_INSTALL_PREFIX})
  ")
elseif(CPACK_GENERATOR MATCHES "TGZ_NODEP")
  set(CPACK_GENERATOR "TGZ")
else(CPACK_GENERATOR MATCHES "DEB")
  message(FATAL_ERROR "Uknown CPack generator")
endif(CPACK_GENERATOR MATCHES "DEB")
