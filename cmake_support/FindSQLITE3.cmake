# Once done this will define
# - Find SQLite3
# Find the native SQLite3 includes and library
# 
#  SQLITE3_INCLUDE_DIRS  - SQLite3 include dir.
#  SQLITE3_LIBRARIES     - List of libraries when using SQLite3.
#  SQLITE3_FOUND         - True if SQLite3 found.

find_path(SQLITE3_INCLUDE_DIR sqlite3.h)
find_library(SQLITE3_LIBRARY NAMES sqlite3)

# handle the QUIETLY and REQUIRED arguments and set SQLITE3_FOUND to TRUE if 
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (SQLITE3 DEFAULT_MSG SQLITE3_LIBRARY SQLITE3_INCLUDE_DIR)

mark_as_advanced(SQLITE3_LIBRARY SQLITE3_INCLUDE_DIR)

SET(SQLITE3_INCLUDE_DIRS ${SQLITE3_INCLUDE_DIR})
SET(SQLITE3_LIBRARIES ${SQLITE3_LIBRARY})
