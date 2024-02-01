# Try to find SpaceNav lib
# adapted from FINDEIGEN.cmake of OpenStructure by Stefan Scheuber <stefan.scheuber@unibas.ch>
# license see FINDEIGEN.cmake 
#
# Once done this will define
#
#  SPNAV_FOUND - system has spacenav lib
#  SPNAV_LIBRARIES - List of libraries when using SpNav.
#  SPNAV_INCLUDE_DIR - the spacenav include directory

if (SPNAV_INCLUDE_DIR)
  # in cache already
  set(SPNAV_FOUND TRUE)

else (SPNAV_INCLUDE_DIR)
  find_path (SPNAV_INCLUDE_DIR NAMES spnav.h)
  find_library (SPNAV_LIBRARIES NAMES spnav)
endif(SPNAV_INCLUDE_DIR)

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (SPNAV DEFAULT_MSG SPNAV_LIBRARIES SPNAV_INCLUDE_DIR)

mark_as_advanced (SPNAV_LIBRARIES SPNAV_INCLUDE_DIR)





