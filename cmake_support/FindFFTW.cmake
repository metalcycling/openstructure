# Try to find FFTW
#
# Once done this will define
# - Find FFTW
# Find the native FFTW includes and library
# 
#  FFTW_INCLUDE_DIRS - FFTW include dirs.
#  FFTW_LIBRARIES    - List of libraries when using FFTW.
#  FFTW_FOUND        - True if FFTW found.

find_path (FFTW_INCLUDE_DIR fftw3.h)
if (_DOUBLE_PREC)
  find_library (FFTW_LIBRARY NAMES fftw3)
  find_library (FFTW_THREADS_LIB NAMES fftw3_threads)
elseif(NOT _DOUBLE_PREC)
  find_library (FFTW_LIBRARY NAMES fftw3f)
  find_library (FFTW_THREADS_LIBRARY NAMES fftw3f_threads)
endif()

if(FFTW_THREADS_LIBRARY)
  set(FFTW_USE_THREADS TRUE)
  set(FFTW_LIBRARIES ${FFTW_LIBRARY} ${FFTW_THREADS_LIBRARY} )
else(FFTW_THREADS_LIBRARY)
  set(FFTW_LIBRARIES ${FFTW_LIBRARY} )
endif(FFTW_THREADS_LIBRARY)

set(FFTW_INCLUDE_DIRS ${FFTW_INCLUDE_DIR})

# handle the QUIETLY and REQUIRED arguments and set FFTW_FOUND to TRUE if 
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (FFTW DEFAULT_MSG FFTW_LIBRARY FFTW_INCLUDE_DIR)

mark_as_advanced (FFTW_LIBRARY FFTW_THREADS_LIBRARY FFTW_INCLUDE_DIR)
