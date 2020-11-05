#-------------------------------------------------------------------------------
# Check for Python Libraries
#
#    PYTHON_IGNORE_FRAMEWORKS if set, do not check for python frameworks. 
#                             has meaning on macOS only
#    PYTHON_ROOT              Prefix for python libraries
#    PYTHON_MIN_VERSION       minimal python version required
#  
# When Python is found, the result is placed in the following variables:
# 
#    PYTHON_LIBRARIES         is set to the library and linker flags used to
#                             link against python
#    PYTHON_VERSION           is set to the version of python
#    PYTHON_INCLUDE_PATH      is set to the path that contains Python.h
#    PYTHON_BINARY            is set to the path to the python executable
#
# Author: Marco Biasini
#-------------------------------------------------------------------------------

set(PYTHON_VERSIONS 3.9 3.8 3.7 3.6 )
set(PYTHON_MIN_VERSION 3.6.0)

#-------------------------------------------------------------------------------
# check for python framework
# this macro honours the values of PYTHON_ROOT
#-------------------------------------------------------------------------------
macro(check_for_python_framework)
  set(_FRAMEWORK_SEARCH_PATHS /Library/Frameworks/ /System/Library/Frameworks)
  if(PYTHON_ROOT)
    set(_FRAMEWORK_SEARCH_PATHS ${PYTHON_ROOT}/Library/Frameworks)
  endif()
  foreach(_PATH ${_FRAMEWORK_SEARCH_PATHS})
    set(_FULL_FRAMEWORK_NAME "${_PATH}/Python.framework")
    if(EXISTS ${_FULL_FRAMEWORK_NAME})
      set(PYTHON_FRAMEWORK ON)
      set(PYTHON_INCLUDE_PATH "${_FULL_FRAMEWORK_NAME}/Headers")
      set(PYTHON_FRAMEWORK_PATH "${_FULL_FRAMEWORK_NAME}/Python")
    endif()
  endforeach()
endmacro()



macro(_find_python PYTHON_ROOT VERSION)
  string(REPLACE "." "" _VERSION_NO_DOTS "${VERSION}")
  if(PYTHON_ROOT)
    find_library(PYTHON_LIBRARIES
      NAMES "python${_VERSION_NO_DOTS}" "python${VERSION}" 
      "python${_VERSION_NO_DOTS}m" "python${VERSION}m"
      HINTS "${PYTHON_ROOT}"
      PATH_SUFFIXES lib libs
      NO_SYSTEM_ENVIRONMENT_PATH NO_DEFAULT_PATH
    )
    find_path(PYTHON_INCLUDE_PATH
      NAMES Python.h
      HINTS "${PYTHON_ROOT}/include"
      PATH_SUFFIXES include "python${_VERSION_NO_DOTS}" "python${VERSION}" 
      "python${_VERSION_NO_DOTS}m" "python${VERSION}m"
      NO_SYSTEM_ENVIRONMENT_PATH NO_DEFAULT_PATH
    )
  else()
    find_library(PYTHON_LIBRARIES
      NAMES "python${_VERSION_NO_DOTS}" "python${VERSION}" 
      "python${_VERSION_NO_DOTS}m" "python${VERSION}m"
      PATH_SUFFIXES lib
    )
    find_path(PYTHON_INCLUDE_PATH
      NAMES Python.h
      PATH_SUFFIXES include "python${_VERSION_NO_DOTS}" "python${VERSION}" 
      "python${_VERSION_NO_DOTS}m" "python${VERSION}m" 
    )    
  endif()  
endmacro()

macro(_find_python_bin PYTHON_ROOT VERSION)
  string(REPLACE "." "" _VERSION_NO_DOTS "${VERSION}")
  if(PYTHON_ROOT)
    find_program(PYTHON_BINARY
      NAMES "python${_VERSION_NO_DOTS}" "python${VERSION}" python.exe
      HINTS "${PYTHON_ROOT}"
      PATH_SUFFIXES bin
      NO_SYSTEM_ENVIRONMENT_PATH NO_DEFAULT_PATH
    )
  else()
    find_program(PYTHON_BINARY
      NAMES "python${_VERSION_NO_DOTS}" "python${VERSION}"
      HINTS "${CMAKE_PREFIX_PATH}"
      PATH_SUFFIXES bin
    )  
  endif()  
endmacro()

#-------------------------------------------------------------------------------
# check for python lib
#
# this macro honours the values of PYTHON_ROOT and PYTHON_VERSION
#-------------------------------------------------------------------------------
macro(check_for_python_lib)
  if(PYTHON_VERSION)
    _find_python("${PYTHON_ROOT}" "${PYTHON_VERSION}")
  else()
    foreach(_VERSION ${PYTHON_VERSIONS})
     if((${PYTHON_MIN_VERSION} VERSION_LESS ${_VERSION}) OR
        (${PYTHON_MIN_VERSION} VERSION_EQUAL ${_VERSION}))
        _find_python("${PYTHON_ROOT}" "${_VERSION}")
        if(PYTHON_LIBRARIES)
          set(PYTHON_VERSION "${_VERSION}")
          break()
        endif()
      endif()
    endforeach()
  endif()
  # fallback to non-versioned naming scheme
  if (NOT $PYTHON_LIBRARIES)
    _find_python("${PYTHON_ROOT}" "")
  endif()
endmacro()

macro(check_for_python_binary)
  if(PYTHON_VERSION)
    _find_python_bin("${PYTHON_ROOT}" "${PYTHON_VERSION}")
  else()
    foreach(_VERSION ${PYTHON_VERSIONS})
      if((${PYTHON_MIN_VERSION} VERSION_LESS ${_VERSION}) OR
         (${PYTHON_MIN_VERSION} VERSION_EQUAL ${_VERSION}))
        _find_python_bin("${PYTHON_ROOT}" "${_VERSION}")
        if(PYTHON_BINARY)
          set(PYTHON_VERSION "${_VERSION}")
          # disallow all versions except for the one we just found. This makes
          # sure we don't mismatch the python binary and the libraries.
          set(PYTHON_VERSIONS "${_VERSION}")
          break()
        endif()
      endif()
    endforeach()
  endif()
  if (NOT PYTHON_BINARY)
    _find_python("${PYTHON_ROOT}" "")
  endif()
endmacro()

if(NOT PYTHON_ROOT)
  if(WIN32)
    set(PYTHON_ROOT "${CMAKE_PREFIX_PATH}")
  else()
    set(PYTHON_ROOT "/usr")
  endif()
endif()
if(APPLE AND NOT PYTHON_IGNORE_FRAMEWORKS)
  check_for_python_framework()
endif()

# first check for python binary.
check_for_python_binary()

if(NOT PYTHON_FRAMEWORK_FOUND)
  check_for_python_lib()
endif()

mark_as_advanced(
  PYTHON_LIBRARIES
  PYTHON_INCLUDE_PATH
  PYTHON_VERSION
  PYTHON_BINARY
)

if(PYTHON_LIBRARIES)
  if(PYTHON_FRAMEWORK)
    set(PYTHON_LIBRARIES "${PYTHON_FRAMEWORK_PATH}"
        CACHE FILEPATH "Python Libraries" FORCE)
  else()
    set(PYTHON_LIBRARIES "${PYTHON_LIBRARIES}"
        CACHE FILEPATH "Python Libraries" FORCE)
  endif()
  set(PYTHON_INCLUDE_PATH "${PYTHON_INCLUDE_PATH}"
      CACHE FILEPATH "Python Include Path" FORCE)
endif()

if (PYTHON_BINARY)
  set(PYTHON_VERSION "${PYTHON_VERSION}"
      CACHE STRING "Python Version" FORCE)
  set(PYTHON_BINARY "${PYTHON_BINARY}"
      CACHE FILEPATH "Python Binary" FORCE)
endif()
