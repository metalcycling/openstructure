# - Try to find Eigen3 lib
# Once done this will define
#
#  EIGEN3_FOUND - system has eigen lib
#  EIGEN3_INCLUDE_DIRS - the eigen include directories

find_path(EIGEN3_INCLUDE_DIR NAMES Eigen/Core
PATH_SUFFIXES eigen3
)

set(EIGEN3_INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR}) 

include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (EIGEN DEFAULT_MSG EIGEN3_INCLUDE_DIR)

mark_as_advanced (EIGEN3_INCLUDE_DIR)





