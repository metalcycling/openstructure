#-------------------------------------------------------------------------------
#  Author:  Andreas Schenk
#
#  This file contains the CPack configuration common to all platforms.
#
#-------------------------------------------------------------------------------

#only include deployment configureation for the main project and not for subprojects
if(${PROJECT_SOURCE_DIR} STREQUAL ${CMAKE_SOURCE_DIR})

  INCLUDE(InstallRequiredSystemLibraries)
  
  set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "OpenStructure Framework")
  set(CPACK_PACKAGE_VENDOR "2008-2020 -- The OpenStructure Team")
  set(CPACK_PACKAGE_DESCRIPTION_FILE "${PROJECT_SOURCE_DIR}/ReadMe.txt")
  set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE.txt")
  set(CPACK_PACKAGE_VERSION_MAJOR ${OST_VERSION_MAJOR})
  set(CPACK_PACKAGE_VERSION_MINOR ${OST_VERSION_MINOR})
  set(CPACK_PACKAGE_VERSION_PATCH ${OST_VERSION_PATCH})
  
  # using include instead of add_subdirectory avoid creating a child scope
  if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    include(macos/CMakeLists.txt) 
  elseif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    include(linux/CMakeLists.txt) 
  elseif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    include(win/CMakeLists.txt) 
  endif()

  INCLUDE(CPack) 

endif(${PROJECT_SOURCE_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
