#-------------------------------------------------------------------------------
#  Author:  Andreas Schenk
#
#  This file contains the Windows specific CPack configuration.
#
#-------------------------------------------------------------------------------

install(CODE "
  if(NOT ${CMAKE_INSTALL_PREFIX} STREQUAL \${CMAKE_INSTALL_PREFIX})
    message(FATAL_ERROR "The Windows installer cannot yet be created with CPack.")            
  endif(NOT ${CMAKE_INSTALL_PREFIX} STREQUAL \${CMAKE_INSTALL_PREFIX})
")

