# Try to find OpenMM
#
# Once done this will define
# - Find OpenMM
# Find the native OpenMM includes and library
# 
#  OPEN_MM_INCLUDE_DIR  - OpenMM include dirs.
#  OPEN_MM_LIBRARIES    - List of libraries when using OpenMM.
#  OPEN_MM_FOUND        - True if OpenMM found.
#  OPEN_MM_PLUGIN_DIR   - Default for ressource plugins

if (OPEN_MM_INCLUDE_DIR)
  set(OPEN_MM_FOUND TRUE)
else (OPEN_MM_INCLUDE_DIR)
  find_path (OPEN_MM_INCLUDE_DIR OpenMM.h)
endif (OPEN_MM_INCLUDE_DIR)

if (OPEN_MM_LIBRARY)
  set(OPEN_MM_FOUND TRUE)
else (OPEN_MM_LIBRARY)
  find_library (OPEN_MM_LIBRARY NAMES OpenMM)
endif(OPEN_MM_LIBRARY)

set(OPEN_MM_LIBRARIES ${OPEN_MM_LIBRARY})
set(OPEN_MM_INCLUDE_DIRS ${OPEN_MM_INCLUDE_DIR})

# handle the QUIETLY and REQUIRED arguments and set OpenMM_FOUND to TRUE if 
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (OpenMM DEFAULT_MSG OPEN_MM_LIBRARY OPEN_MM_INCLUDE_DIR)

if (OPEN_MM_LIBRARY)
  set(OPEN_MM_FOUND TRUE)
  if (NOT OPEN_MM_PLUGIN_DIR)
    set(OPEN_MM_PLUGIN_DIR "/usr/local/openmm/lib/plugins")
  endif (NOT OPEN_MM_PLUGIN_DIR)
  if (NOT EXISTS "${OPEN_MM_PLUGIN_DIR}")
    message(FATAL_ERROR "OpenMM plugin directory does not exist: '${OPEN_MM_PLUGIN_DIR}'")
  endif (NOT EXISTS "${OPEN_MM_PLUGIN_DIR}")
endif (OPEN_MM_LIBRARY)

mark_as_advanced (OPEN_MM_LIBRARY OPEN_MM_INCLUDE_DIR OPEN_MM_FOUND
                  OPEN_MM_PLUGIN_DIR)
