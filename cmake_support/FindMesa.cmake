set(OPENGL_FOUND FALSE)

if (OPENGL_INCLUDE_DIR)

  # in cache already
  set(OPENGL_FOUND TRUE)

else (OPENGL_INCLUDE_DIR)

  if(DEFINED MESA_PREFIX)
    set(_MESA_PREFIX ${MESA_PREFIX})
  elseif(DEFINED ENV{MESA_PREFIX})
    set(_MESA_PREFIX $ENV{MESA_PREFIX})
  else()
    message(FATAL_ERROR "no MESA_PREFIX given or found in environment")
  endif()

  if(DEFINED _MESA_PREFIX)
    set(OPENGL_INCLUDE_DIR "${_MESA_PREFIX}/include" CACHE PATH "Mesa include")
    set(OPENGL_LIBRARY_DIR "${_MESA_PREFIX}/lib" CACHE PATH "Mesa lib")
    set(OPENGL_gl_LIBRARY "${_MESA_PREFIX}/lib/libGL.so" CACHE PATH "Mesa gl lib")
    set(OPENGL_glu_LIBRARY "${_MESA_PREFIX}/lib/libGLU.so" CACHE PATH "Mesa glu lib")

    set(OPENGL_LIBRARIES ${_MESA_PREFIX}/lib/libGL.so ${_MESA_PREFIX}/lib/libGLU.so ${_MESA_PREFIX}/lib/libOSMesa.so CACHE PATH "Mesa libs")
    set(OPENGL_FOUND TRUE)
    message(STATUS "Found Mesa include dir: ${OPENGL_INCLUDE_DIR}")
    message(STATUS "Found Mesa lib dir: ${OPENGL_LIBRARY_DIR}")
    message(STATUS "Using Mesa libs: ${OPENGL_LIBRARIES}")
  endif()

endif()

if(NOT OPENGL_FOUND)
    message(FATAL_ERROR "no Mesa found, please provide the correct location with MESA_PREFIX")
endif()

mark_as_advanced (OPENGL_INCLUDE_DIR)
