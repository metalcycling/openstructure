#-------------------------------------------------------------------------------
# Authors: Marco Biasini, Juergen Haas, Andreas Schenk
#
# This file contains a bunch of useful macros to facilitate the build-system
# configuration for the modules.
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# map macro
#
# this function emulates a map/dict data type
#-------------------------------------------------------------------------------

function(map COMMAND MAPNAME)
  set (_KEYS ${MAPNAME}_MAP_KEYS )
  set (_VALUES ${MAPNAME}_MAP_VALUES)
  if(${COMMAND} STREQUAL SET)
    list(REMOVE_AT ARGN 0)
    list(FIND ${_KEYS} ${ARGV2} _MAP_INDEX)
    if(_MAP_INDEX EQUAL -1)
      list(APPEND ${_KEYS} ${ARGV2})
      set(${_KEYS} ${${_KEYS}} PARENT_SCOPE)
      set(${_VALUES}_${ARGV2}  ${ARGN} PARENT_SCOPE)
    else()
      set(${_VALUES}_${ARGV2}  ${ARGN} PARENT_SCOPE)
    endif()
  elseif(${COMMAND} STREQUAL GET)
    list(FIND ${_KEYS} ${ARGV2} _MAP_INDEX)
    if(_MAP_INDEX EQUAL -1)
      MESSAGE(FATAL_ERROR "Unknown key: " ${ARGV2})
    endif()
    set(${ARGV3} ${${_VALUES}_${ARGV2}} PARENT_SCOPE)
  elseif(${COMMAND} STREQUAL KEYS)
    set(${ARGV2} ${${_KEYS}} PARENT_SCOPE)
  elseif(${COMMAND} STREQUAL CREATE)
    set(${_KEYS}  "" PARENT_SCOPE)
  elseif(${COMMAND} STREQUAL LENGTH)
    list(LENGTH ${_KEYS} _L)
    set(${ARGV2} ${_L} PARENT_SCOPE)
  else()
    MESSAGE(FATAL_ERROR "Unknown map command:" ${COMMAND})
  endif()
endfunction()


#-------------------------------------------------------------------------------
# check_architecture
#
# detect architecture based on void pointer size. the output is stored in the
# 3 variables OS_32_BITS, OS_64_BITS and CMAKE_NATIVE_ARCH. The former two are
# set to 0 and 1 accordingly and CMAKE_NATIVE_ARCH is set to the 32 or 64 bit.
#-------------------------------------------------------------------------------
macro(check_architecture)
  include(CheckTypeSize)
  check_type_size(void*  SIZEOF_VOID_PTR)
  if(${SIZEOF_VOID_PTR} MATCHES "^8$")
    set(OS_32_BITS 0)
    set(OS_64_BITS 1)
    set(CMAKE_NATIVE_ARCH 64)
  else()
    set(OS_32_BITS 1)
    set(OS_64_BITS 0)
    set(CMAKE_NATIVE_ARCH 32)
  endif()
endmacro()

#-------------------------------------------------------------------------------
# this macro has been adapted from
# http://www.cmake.org/Wiki/CMakeMacroParseArguments
#-------------------------------------------------------------------------------
macro(parse_argument_list PREFIX ARG_NAMES OPT_NAMES)
  set(_DEFAULT_ARGS)
  # reset variables
  foreach(${_ARG_NAME} ${ARG_NAMES})
    set(PREFIX_${_ARG_NAME})
  endforeach()
  foreach(_OPT_NAME ${OPT_NAMES})
    set(PREFIX_${_OPT_NAME} FALSE)
  endforeach()
  set(_CURR_ARG_NAME DEF_ARGS)
  set(_CURR_ARG_LIST)
  # loop over parameter list and split by ARG_NAMES
  foreach(_ARG ${ARGN})
    set(_LARG_NAMES ${ARG_NAMES})  
    list(FIND _LARG_NAMES ${_ARG} _IS_ARG_NAME)
    if (_IS_ARG_NAME GREATER -1)
      set(${PREFIX}_${_CURR_ARG_NAME} ${_CURR_ARG_LIST})
      set(_CURR_ARG_NAME "${_ARG}")
      set(_CURR_ARG_LIST)
    else()
    set(_LOPT_NAMES ${OPT_NAMES})  
      list(FIND _LOPT_NAMES ${_ARG} _IS_OPT_NAME)
      if (_IS_OPT_NAME GREATER -1)
        set(${PREFIX}_${_ARG} TRUE)
      else()
        list(APPEND _CURR_ARG_LIST "${_ARG}")
      endif()
    endif()
  endforeach()
  set(${PREFIX}_${_CURR_ARG_NAME} ${_CURR_ARG_LIST})
endmacro()

#-------------------------------------------------------------------------------
# copy_if_different
#
# copies file from source directory to destination directory, but only if its
# content changed.
#-------------------------------------------------------------------------------
macro(copy_if_different FROM_DIR TO_DIR FILES TARGETS TARGET)
  foreach(SRC ${FILES})
      set(SRCFILE ${SRC})
      if("${FROM_DIR}" STREQUAL "" OR "${FROM_DIR}" STREQUAL "./")
          set(FROM ${SRC})
      else()
          set(FROM ${FROM_DIR}/${SRC})
      endif()
      if("${TO_DIR}" STREQUAL "")
          set(TO ${SRCFILE})
      else()
          get_filename_component(TOFILE ${SRC} NAME)      
          set(TO ${TO_DIR}/${TOFILE})
      endif()
      file(MAKE_DIRECTORY  ${TO_DIR})
      add_custom_command(TARGET "${TARGET}" PRE_BUILD
          DEPENDS ${FROM}
          COMMAND ${CMAKE_COMMAND} -E copy_if_different ${FROM} ${TO})
  endforeach()
endmacro()


#-------------------------------------------------------------------------------
# parse_file_list
#
# this macro splits a list of files with IN_DIR statements and fills them into a map
# where the key is the directory name
#-------------------------------------------------------------------------------
macro(parse_file_list FILELIST FILEMAP)
  set(_EXPECT_IN_DIR FALSE)
  map(CREATE ${FILEMAP})
  set(_CURRENT_LIST)
  foreach(_ITEM ${FILELIST})
    if (_ITEM STREQUAL "IN_DIR")
      set(_EXPECT_IN_DIR TRUE)
    else()
      if (_EXPECT_IN_DIR)
        set(_EXPECT_IN_DIR FALSE)
        map(SET ${FILEMAP} ${_ITEM} ${_CURRENT_LIST})
        set(_CURRENT_LIST)
      else()
        list(APPEND _CURRENT_LIST "${_ITEM}")
      endif()
    endif()
  endforeach()
  if(_CURRENT_LIST)
    map(SET ${FILEMAP} "." ${_CURRENT_LIST})
  endif()
endmacro()


#-------------------------------------------------------------------------------
# Synopsis:
#   module(NAME name SOURCES source1 source2 HEADERS header1 header2 
#          [IN_DIR dir] [header3 header4 [IN_DIR dir]] [DEPENDS_ON dep1 dep2]
#          [HEADER_OUTPUT_DIR dir]
#          [LINK link_cmds])
# Description:
#   Define an OpenStructure module.
#-------------------------------------------------------------------------------
macro(module)
  #-----------------------------------------------------------------------------
  # deal with arguments
  #-----------------------------------------------------------------------------
  set(_ARGS "NAME;SOURCES;HEADERS;DEPENDS_ON;LINK;HEADER_OUTPUT_DIR;PREFIX")
  set(_ARG_PREFIX ost)  
  parse_argument_list(_ARG "${_ARGS}" "NO_STATIC" ${ARGN})  
  if (NOT _ARG_NAME)
    message(FATAL_ERROR 
            "invalid use of module(): a module name must be provided")
  endif()

  if (ENABLE_STATIC AND _ARG_NO_STATIC)
    return()
  endif()
  if (_ARG_HEADER_OUTPUT_DIR)
    set(_HEADER_OUTPUT_DIR ${_ARG_HEADER_OUTPUT_DIR})
  else()
    set(_HEADER_OUTPUT_DIR "${_ARG_PREFIX}/${_ARG_NAME}")
  endif()
  set(_LIB_NAME ${_ARG_PREFIX}_${_ARG_NAME})
  string(TOUPPER ${_LIB_NAME} _UPPER_LIB_NAME)  
  #-----------------------------------------------------------------------------
  # create library  
  #-----------------------------------------------------------------------------
  file(MAKE_DIRECTORY ${LIB_STAGE_PATH})
  file(MAKE_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
  file(MAKE_DIRECTORY ${LIBEXEC_STAGE_PATH})
  file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/tests")
  if (NOT TARGET create_stage)
    add_custom_target(create_stage COMMAND ${CMAKE_COMMAND} -E make_directory ${LIB_STAGE_PATH}
                                   COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
                                   COMMAND ${CMAKE_COMMAND} -E make_directory ${LIBEXEC_STAGE_PATH}
                                   COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_BINARY_DIR}/tests")
  endif()
  if (_ARG_SOURCES)
    # when there is at least one source file, we build a library
    set(_ABS_SOURCE_NAMES)
    foreach(_SOURCE ${_ARG_SOURCES})
      if (IS_ABSOLUTE ${_SOURCE})
        list(APPEND _ABS_SOURCE_NAMES "${_SOURCE}")
      else()
        list(APPEND _ABS_SOURCE_NAMES "${CMAKE_CURRENT_SOURCE_DIR}/${_SOURCE}")
      endif()
    endforeach()
    if (ENABLE_STATIC AND NOT _ARG_NO_STATIC)
      add_library(${_LIB_NAME} STATIC ${_ABS_SOURCE_NAMES})
    else()
      add_library(${_LIB_NAME} SHARED ${_ABS_SOURCE_NAMES})
    endif()
    set_target_properties(${_LIB_NAME} 
                          PROPERTIES OUTPUT_NAME ${_LIB_NAME}
                                     PROJECT_LABEL ${_ARG_NAME}
                                     EchoString   ${_ARG_NAME}
                                     MODULE_DEPS "${_ARG_DEPENDS_ON}")
    get_target_property(_DEFS ${_LIB_NAME} COMPILE_DEFINITIONS)
    add_dependencies(${_LIB_NAME} create_stage)
    set_target_properties(${_LIB_NAME} PROPERTIES
                          COMPILE_DEFINITIONS OST_MODULE_${_UPPER_LIB_NAME})
    set_target_properties(${_LIB_NAME} PROPERTIES
                          VERSION ${OST_VERSION_STRING}
                          SOVERSION ${OST_VERSION_MAJOR}.${OST_VERSION_MINOR})    
    set_target_properties(${_LIB_NAME} PROPERTIES
                          LIBRARY_OUTPUT_DIRECTORY ${LIB_STAGE_PATH}
                          ARCHIVE_OUTPUT_DIRECTORY ${LIB_STAGE_PATH}
                          RUNTIME_OUTPUT_DIRECTORY ${LIB_STAGE_PATH})
    if (APPLE)
      set_target_properties(${_LIB_NAME} PROPERTIES
                            LINK_FLAGS "-Wl,-rpath,@loader_path"
                            INSTALL_NAME_DIR "@rpath")
    endif()
    if (WIN32)
      install(TARGETS ${_LIB_NAME} RUNTIME DESTINATION bin)
    else()
      if (ENABLE_STATIC)
        install(TARGETS ${_LIB_NAME} ARCHIVE DESTINATION "${LIB_DIR}")
      else()
        install(TARGETS ${_LIB_NAME} LIBRARY DESTINATION "${LIB_DIR}")
      endif()
    endif()                          
    if (_ARG_LINK)
      target_link_libraries(${_LIB_NAME} ${_ARG_LINK})
    endif()
    foreach(_DEPENDENCY ${_ARG_DEPENDS_ON})
      target_link_libraries(${_LIB_NAME} ${_DEPENDENCY})
    endforeach()
    if (ENABLE_STATIC)
      target_link_libraries(${_LIB_NAME} ${STATIC_LIBRARIES})
    endif()
  else()
    add_custom_target("${_LIB_NAME}" ALL)
    set_target_properties("${_LIB_NAME}" PROPERTIES HEADER_ONLY 1 
                          MODULE_DEPS "${_ARG_DEPENDS_ON}")
  endif()
  #-----------------------------------------------------------------------------
  # stage headers  
  #-----------------------------------------------------------------------------
  if (_ARG_HEADERS)
    stage_and_install_headers("${_ARG_HEADERS}" "${_HEADER_OUTPUT_DIR}" "${_LIB_NAME}")
  endif()
endmacro()

#-------------------------------------------------------------------------------
# macro stage_and_install_headers
#-------------------------------------------------------------------------------
macro(stage_and_install_headers HEADERLIST HEADER_OUTPUT_DIR TARGET)
  add_custom_target("${TARGET}_headers" COMMENT "")
  add_dependencies("${TARGET}" "${TARGET}_headers")
  add_dependencies("${TARGET}_headers" create_stage)
  parse_file_list("${HEADERLIST}" _HEADER_MAP)
  map(KEYS _HEADER_MAP _HEADER_MAP_KEYS)
  foreach(_DIR ${_HEADER_MAP_KEYS})
    map(GET _HEADER_MAP ${_DIR} _HEADERS)
    set(_HDR_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${_DIR}")
    set(_ABS_HEADER_NAMES)
    foreach(_HDR ${_HEADERS})
      list(APPEND _ABS_HEADER_NAMES ${_HDR_SOURCE_DIR}/${_HDR})
    endforeach()
    set(_HDR_STAGE_DIR "${HEADER_OUTPUT_DIR}/${_DIR}")
    set(_FULL_HEADER_DIR "${HEADER_STAGE_PATH}/${_HDR_STAGE_DIR}")
    copy_if_different("" "${_FULL_HEADER_DIR}" "${_ABS_HEADER_NAMES}" "" "${TARGET}_headers")
    if(NOT WIN32)
      install(FILES ${_ABS_HEADER_NAMES} DESTINATION "include/${_HDR_STAGE_DIR}")
    endif(NOT WIN32)
  endforeach()
endmacro()


#-------------------------------------------------------------------------------
# Synopsis
#   executable(NAME exe_name SOURCES source1 source2 LINK link1 link2)
#
# Description:
#  Compile, link and stage a C++ executable
#-------------------------------------------------------------------------------
macro(executable)
  parse_argument_list(_ARG 
                      "NAME;SOURCES;LINK;DEPENDS_ON" "NO_RPATH;STATIC" ${ARGN})
  if (NOT _ARG_NAME)
    message(FATAL_ERROR "invalid use of executable(): a name must be provided")
  endif()
  add_executable(${_ARG_NAME} ${_ARG_SOURCES})
  if (APPLE AND NOT _ARG_NO_RPATH AND NOT ENABLE_STATIC)
    set_target_properties(${_ARG_NAME} PROPERTIES
                          LINK_FLAGS "-Wl,-rpath,@loader_path/../lib/")
  endif()
  if (_ARG_LINK)
    target_link_libraries(${_ARG_NAME} ${_ARG_LINK})
  endif()
  foreach(_DEP ${_ARG_DEPENDS_ON})
    target_link_libraries(${_ARG_NAME} ${_DEP})
  endforeach()
  if (ENABLE_STATIC AND _ARG_STATIC)
    target_link_libraries(${_ARG_NAME} ${STATIC_LIBRARIES})
    if (UNIX AND NOT APPLE)
      set_target_properties(${_ARG_NAME} PROPERTIES LINK_SEARCH_START_STATIC TRUE)
      set_target_properties(${_ARG_NAME} PROPERTIES LINK_SEARCH_END_STATIC TRUE)
      if (OST_GCC_LESS_45)
        set_target_properties(${_ARG_NAME} PROPERTIES LINK_FLAGS
                              "-static-libgcc -static -pthread")
      else()
        set_target_properties(${_ARG_NAME} PROPERTIES LINK_FLAGS
                              "-static-libgcc -static-libstdc++ -static -pthread")
      endif()
    endif()
  endif()
  install(TARGETS ${_ARG_NAME} DESTINATION bin)
endmacro()


#-------------------------------------------------------------------------------
# Synopsis
#   executable_libexec(NAME exe_name SOURCES source1 source2 LINK link1 link2)
#
# Description:
#  Compile, link and stage a C++ executable into the libexec directory
#-------------------------------------------------------------------------------
macro(executable_libexec)
  parse_argument_list(_ARG 
                      "NAME;SOURCES;LINK;DEPENDS_ON" "NO_RPATH;STATIC" ${ARGN})
  if (NOT _ARG_NAME)
    message(FATAL_ERROR "invalid use of executable(): a name must be provided")
  endif()   
  add_executable(${_ARG_NAME} ${_ARG_SOURCES})
  set_target_properties(${_ARG_NAME}
                        PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                       "${LIBEXEC_STAGE_PATH}")  
  set_target_properties(${_ARG_NAME}
                        PROPERTIES RUNTIME_OUTPUT_DIRECTORY_RELEASE
                       "${LIBEXEC_STAGE_PATH}")  
  set_target_properties(${_ARG_NAME}
                        PROPERTIES RUNTIME_OUTPUT_DIRECTORY_DEBUG
                       "${LIBEXEC_STAGE_PATH}")  
  if (NOT _ARG_NO_RPATH AND NOT _ARG_STATIC)
    if (APPLE)
      set_target_properties(${_ARG_NAME} PROPERTIES
                            LINK_FLAGS "-Wl,-rpath,@loader_path/../../lib")
    elseif (UNIX)
      set_target_properties(${_ARG_NAME} PROPERTIES INSTALL_RPATH "$ORIGIN/../../${LIB_DIR}")
    endif (APPLE)
  endif (NOT _ARG_NO_RPATH AND NOT _ARG_STATIC)
  if (_ARG_LINK)
    target_link_libraries(${_ARG_NAME} ${_ARG_LINK})
  endif()
  if (ENABLE_STATIC AND _ARG_STATIC)
    target_link_libraries(${_ARG_NAME} ${STATIC_LIBRARIES})
    set_target_properties(${_ARG_NAME}
                          PROPERTIES LINK_SEARCH_END_STATIC TRUE)  

  endif()
  foreach(_DEP ${_ARG_DEPENDS_ON})
    target_link_libraries(${_ARG_NAME} ${_DEP})
  endforeach()
  install(TARGETS ${_ARG_NAME} DESTINATION ${LIBEXEC_PATH})
endmacro()

#-------------------------------------------------------------------------------
# Synopsis:
#   substitute(IN_FILE in_file OUT_FILE out_file DICT a=b c=d)
#
#-------------------------------------------------------------------------------
macro(substitute)
 parse_argument_list(_ARG 
                     "IN_FILE;OUT_FILE;DICT" "" ${ARGN})
  if (NOT _ARG_IN_FILE)
    message(FATAL_ERROR "invalid use of substitute(): no IN_FILE given")
  endif()
  if (NOT _ARG_OUT_FILE)
    message(FATAL_ERROR "invalid use of substitute(): no OUT_FILE given")
  endif()
  set(_SUBST_DICT -DINPUT_FILE=${_ARG_IN_FILE} -DOUT_FILE=${_ARG_OUT_FILE})
  foreach(_SUBST ${_ARG_DICT})
    list(APPEND _SUBST_DICT -D${_SUBST})
  endforeach()
  add_custom_target(subst_${_ARG_OUT_FILE} ALL COMMAND 
                    ${CMAKE_COMMAND} ${_SUBST_DICT}
                    -P ${CMAKE_SOURCE_DIR}/cmake_support/substitute.cmake)
endmacro()

#-------------------------------------------------------------------------------
# Synopsis:
#   script(NAME script_name INPUT input_name SUBSTITUTE key=val key=val
#          [OUTPUT_DIR dir] [TARGET target])
#-------------------------------------------------------------------------------
macro(script)
  set(_ARG_OUTPUT_DIR bin)
  parse_argument_list(_ARG 
                      "NAME;INPUT;SUBSTITUTE;TARGET;OUTPUT_DIR" "" ${ARGN})
  if (NOT _ARG_NAME)
    message(FATAL_ERROR "invalid use of executable(): a name must be provided")
  endif()
  set(_INPUT ${_ARG_NAME})
  if (_ARG_INPUT)
    set(_INPUT ${_ARG_INPUT})
  endif()
  if (_ARG_SUBSTITUTE)
    if (NOT _ARG_INPUT)
      message(FATAL_ERROR "script() can only substitute when INPUT is present.")    
    endif()

    substitute(IN_FILE "${CMAKE_CURRENT_SOURCE_DIR}/${_INPUT}" OUT_FILE ${_ARG_NAME} 
               DICT ${_ARG_SUBSTITUTE})
  endif()
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${_ARG_NAME} DESTINATION ${_ARG_OUTPUT_DIR} 
          PERMISSIONS WORLD_EXECUTE GROUP_EXECUTE OWNER_EXECUTE 
                      WORLD_READ GROUP_READ OWNER_READ)
  copy_if_different("./" "${STAGE_DIR}/${_ARG_OUTPUT_DIR}" 
                    "${_ARG_NAME}" "TARGETS" ${_ARG_TARGET})
  add_dependencies(${_ARG_TARGET} subst_${_ARG_NAME})
endmacro()

#-------------------------------------------------------------------------------
# Synopsis:
#   ui_to_python(libname stagedir[input_file1 ...])
# Description:
#   Calls pyuic on every input file. The resulting python files are stored in
#   the variable with name out_files. Use with care. Has originally been 
#   implemented for pyuic4 and ported to pyuic5 without testing. 
#-------------------------------------------------------------------------------
macro(ui_to_python LIBNAME PYMODDIR STAGEDIR)
  set(_input_files ${ARGN})
  add_custom_target("${LIBNAME}_ui" ALL)
  add_dependencies("_${LIBNAME}" "${LIBNAME}_ui")
  find_program(_PYUIC_EXECUTABLE
    NAMES pyuic5-${PYTHON_VERSION} pyuic5 pyuic pyuic5.bat
    PATHS  ENV PATH 
  )  
  if(NOT _PYUIC_EXECUTABLE)
    message(FATAL_ERROR "Could not find pyuic command in " ${QT_BINARY_DIR} " for python version " ${PYTHON_VERSION})
  endif(NOT _PYUIC_EXECUTABLE)
  set(out_files)
  foreach(input_file ${_input_files})
    get_filename_component(_out_file ${input_file} NAME_WE)
    get_filename_component(_in_file ${input_file} ABSOLUTE)
    set(_out_file ${_out_file}_ui.py)
    set(_abs_out_file ${STAGEDIR}/${_out_file})
    add_custom_command(TARGET ${LIBNAME}_ui
                       COMMAND ${_PYUIC_EXECUTABLE} -o ${_abs_out_file} ${_in_file}
                       VERBATIM DEPENDS ${input_file}
                       )
    list(APPEND out_files ${_abs_out_file})
  endforeach()
  compile_py_files(_${LIBNAME} ${STAGEDIR} compiled_files ${out_files})
  install(FILES ${out_files} DESTINATION "${LIB_DIR}/${PYMODDIR}")
  install(FILES ${compiled_files} DESTINATION "${LIB_DIR}/${PYMODDIR}")
endmacro()

#-------------------------------------------------------------------------------
# Synopsis:
#   compile_py_files(module out_dir compiled_files [input_file1 ...])
# Description:
#   Calls pyuic on every input file. The resulting python files are stored in
#   the variable with name compiled_files.
#-------------------------------------------------------------------------------
macro(compile_py_files module out_dir compiled_files_name)
  set(_input_files ${ARGN})
  set(${compiled_files_name})
  foreach(input_file ${_input_files})
    get_filename_component(_out_file ${input_file} NAME_WE)
    get_filename_component(_in_file ${input_file} ABSOLUTE)
    set(_out_file ${out_dir}/${_out_file}.pyc)
    list(APPEND ${compiled_files_name} ${_out_file})
    get_filename_component(_in_name ${input_file} NAME)
    file(MAKE_DIRECTORY  ${out_dir})
    add_custom_command(TARGET ${module}
                       COMMAND ${PYTHON_BINARY} -c "import py_compile;py_compile.compile(\"${_in_file}\",\"${_out_file}\",\"${_in_name}\",doraise=True)"
                       VERBATIM DEPENDS ${input_file}
                       )
  endforeach()
endmacro()

#-------------------------------------------------------------------------------
# Synopsis:
#   pymod(NAME name CPP source1 source2 PY source source2 [IN_DIR dir] 
#         source3 source4 [IN_DIR dir] [LINK link] [OUTPUT_DIR dir] [UI user_interface_files])
#
# Description:
#  Define a python module consisting of C++ type wrappers and/or code written in 
#  Python. user_interface_files are handled with ui_to_python that has blindly
#  been ported from QT4 to QT5 without testing. Use with care.
# NAME is the name of
#-------------------------------------------------------------------------------
macro(pymod)
  #-----------------------------------------------------------------------------
  # deal with arguments
  #-----------------------------------------------------------------------------
  set(_ARG_PREFIX ost)
  parse_argument_list(_ARG 
                      "NAME;CPP;PY;LINK;OUTPUT_DIR;UI;PREFIX" "" ${ARGN})
  if (NOT _ARG_NAME)
    message(FATAL_ERROR "invalid use of pymod(): a name must be provided")
  endif()
  if (ENABLE_STATIC)
    return()
  endif()
  if (_ARG_OUTPUT_DIR)
    set(PYMOD_DIR "python${PYTHON_VERSION}/site-packages/${_ARG_OUTPUT_DIR}")
  else()
    set(PYMOD_DIR "python${PYTHON_VERSION}/site-packages/${_ARG_PREFIX}/${_ARG_NAME}")
  endif()
  set(_LIB_NAME ${_ARG_PREFIX}_${_ARG_NAME})
  set(PYMOD_STAGE_DIR "${LIB_STAGE_PATH}/${PYMOD_DIR}")
  file(MAKE_DIRECTORY ${PYMOD_STAGE_DIR})
  include_directories(${PYTHON_INCLUDE_PATH})
  #-----------------------------------------------------------------------------
  # compile and link C++ wrappers
  #-----------------------------------------------------------------------------
  if (_ARG_CPP)
    add_library("_${_LIB_NAME}" MODULE ${_ARG_CPP})
    set_target_properties("_${_LIB_NAME}"
                          PROPERTIES ECHO_STRING
                          "Building Python Module ${_ARG_NAME}")
    if (_ARG_PREFIX)
      set(_PARENT_NAME "${_ARG_PREFIX}_${_ARG_NAME}")
    else()
      set(_PARENT_NAME "${_ARG_NAME}")
    endif()
    get_target_property(_CUSTOM_CHECK "${_PARENT_NAME}" HEADER_ONLY)
    if (NOT _CUSTOM_CHECK)
      set(_PARENT_LIB_NAME "${_PARENT_NAME}")
    endif()
    target_link_libraries("_${_LIB_NAME}" ${_PARENT_LIB_NAME} 
                          ${PYTHON_LIBRARIES} ${BOOST_PYTHON_LIBRARIES})

    set_target_properties("_${_LIB_NAME}"
                          PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PYMOD_STAGE_DIR})
    set_target_properties("_${_LIB_NAME}"
                          PROPERTIES LIBRARY_OUTPUT_DIRECTORY_RELEASE ${PYMOD_STAGE_DIR})
    set_target_properties("_${_LIB_NAME}"
                          PROPERTIES LIBRARY_OUTPUT_DIRECTORY_DEBUG ${PYMOD_STAGE_DIR})

    if (NOT ENABLE_STATIC)
      if (_USE_RPATH)
        string(REGEX REPLACE "/[^/]*" "/.." inv_pymod_path "/${PYMOD_DIR}")
        set_target_properties("_${_LIB_NAME}"
                              PROPERTIES INSTALL_RPATH "$ORIGIN${inv_pymod_path}/")
      else()
        set_target_properties("_${_LIB_NAME}"
                              PROPERTIES INSTALL_RPATH "")
      endif()
    endif()
    if (APPLE)
      file(RELATIVE_PATH _REL_PATH "${PYMOD_STAGE_DIR}" "${LIB_STAGE_PATH}")
      set_target_properties("_${_LIB_NAME}" PROPERTIES
                            LINK_FLAGS "-Wl,-rpath,@loader_path/${_REL_PATH}"
                            INSTALL_NAME_DIR "@rpath")
    endif()                          
    if (NOT WIN32)
      set_target_properties("_${_LIB_NAME}"
                          PROPERTIES PREFIX "")
    else ()
      set_target_properties("_${_LIB_NAME}"
                          PROPERTIES SUFFIX ".pyd")
    endif()
    install(TARGETS "_${_LIB_NAME}" LIBRARY DESTINATION
            "${LIB_DIR}/${PYMOD_DIR}")
  else()
    add_custom_target("_${_LIB_NAME}" ALL)
  endif()
  #-----------------------------------------------------------------------------
  # build ui files
  #-----------------------------------------------------------------------------
  if (_ARG_UI)
    ui_to_python(${_LIB_NAME} ${PYMOD_DIR} ${PYMOD_STAGE_DIR} ${_ARG_UI})
  endif()  
  #-----------------------------------------------------------------------------
  # compile python files
  #-----------------------------------------------------------------------------
  if (_ARG_PY)
    parse_file_list("${_ARG_PY}" _PYFILE_MAP)
    map(KEYS _PYFILE_MAP _PYFILE_MAP_KEYS)
    foreach(_DIR ${_PYFILE_MAP_KEYS})
      map(GET _PYFILE_MAP ${_DIR} _PY_FILES)
      set(_ABS_PY_FILES)
      set(_PY_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${_DIR}")
      foreach(_PY ${_PY_FILES})
        list(APPEND _ABS_PY_FILES "${_PY_SOURCE_DIR}/${_PY}")
      endforeach()
      install(FILES ${_ABS_PY_FILES} DESTINATION "${LIB_DIR}/${PYMOD_DIR}/${_DIR}")
      string(REPLACE "/" "_" _DIR_NO_SLASH "${_DIR}")
      set(_PYMOD_TARGET "${_LIB_NAME}_${_DIR_NO_SLASH}_pymod")
      string(REPLACE "_." "" _PYMOD_TARGET "${_PYMOD_TARGET}")
      add_custom_target(${_PYMOD_TARGET} ALL)
      copy_if_different("./" "${PYMOD_STAGE_DIR}/${_DIR}"
                        "${_ABS_PY_FILES}" "TARGETS"
                        "${_PYMOD_TARGET}")
      compile_py_files(_${_LIB_NAME} ${PYMOD_STAGE_DIR}/${_DIR} compiled_files ${_ABS_PY_FILES})
      install(FILES ${compiled_files} DESTINATION "${LIB_DIR}/${PYMOD_DIR}/${_DIR}")
    endforeach()
  endif()  
  get_target_property(_MOD_DEPS "${_PARENT_NAME}" MODULE_DEPS)
  if(_MOD_DEPS)
    foreach(dep ${_MOD_DEPS})
       add_dependencies("_${_LIB_NAME}" "_${dep}")
    endforeach()
  endif()
endmacro()

add_custom_target(check)
add_custom_target(check_xml)
if (WIN32)
  set_target_properties(check PROPERTIES EXCLUDE_FROM_ALL "1")
  set_target_properties(check_xml PROPERTIES EXCLUDE_FROM_ALL "1")
endif()

#-------------------------------------------------------------------------------
# ost_unittest
#
# define a unit test
#-------------------------------------------------------------------------------
macro(ost_unittest)
  set(_ARG_PREFIX ost)
  parse_argument_list(_ARG 
                      "MODULE;PREFIX;SOURCES;LINK" "" ${ARGN})
    set(_SOURCES ${_ARG_SOURCES})
    set(CPP_TESTS)
    set(PY_TESTS)
    set(CMAKE_CURRENT_BINARY_DIR "${CMAKE_BINARY_DIR}/tests")
    foreach(src ${_SOURCES})
      if (${src} MATCHES "\\.py$")
       list(APPEND PY_TESTS "${src}")      
      else()
        list(APPEND CPP_TESTS "${CMAKE_CURRENT_SOURCE_DIR}/${src}")
     endif()
    endforeach()
    set(_SOURCES ${CPP_TESTS})
    set(_test_name "${_ARG_PREFIX}_${_ARG_MODULE}_tests")
    if(DEFINED CPP_TESTS)
      if(COMPILE_TESTS)
        add_executable(${_test_name} ${_SOURCES})
      else()
        add_executable(${_test_name} EXCLUDE_FROM_ALL ${_SOURCES})
      endif()
      set_target_properties(${_test_name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests"  )
      set_target_properties(${_test_name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}/tests"  )
      set_target_properties(${_test_name} PROPERTIES RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/tests"  )

      target_link_libraries(${_test_name} ${BOOST_UNIT_TEST_LIBRARIES} "${_ARG_PREFIX}_${_ARG_MODULE}")
      if (WIN32)
        set(TEST_COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${_test_name}.exe || echo)
        message(${TEST_COMMAND})
        set(ENV{PP} ${STAGE_DIR}/bin)
        add_custom_target("${_test_name}_run" 
                        COMMAND echo %PP%
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        COMMENT "running checks for module ${_ARG_MODULE}"
                        DEPENDS ${_test_name})
        add_test("${_test_name}" ${CMAKE_CURRENT_BINARY_DIR}/${_test_name}.exe)
        set_target_properties("${_test_name}_run" PROPERTIES EXCLUDE_FROM_ALL "1")
        add_dependencies(check "${_test_name}_run")
      else()
        add_custom_target("${_test_name}_run"
                        COMMAND OST_ROOT=${STAGE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/${_test_name}
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        COMMENT "running checks for module ${_ARG_MODULE}"
                        DEPENDS ${_test_name})
        add_custom_target("${_test_name}_run_xml"
                        COMMAND OST_ROOT=${STAGE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/${_test_name} --log_format=xml --log_level=all > ${_test_name}_log.xml || echo
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        COMMENT "running checks for module ${_ARG_MODULE}"
                        DEPENDS ${_test_name})
        add_test("${_test_name}" ${CMAKE_CURRENT_BINARY_DIR}/${_test_name} )
        add_dependencies(check_xml "${_test_name}_run_xml")
        add_dependencies(check "${_test_name}_run")
      endif()
      
      if (_ARG_LINK)
        target_link_libraries("${_test_name}" ${_ARG_LINK})
      endif()

      set_target_properties(${_test_name}
                            PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                            "${CMAKE_CURRENT_BINARY_DIR}")
    endif()

    foreach(py_test ${PY_TESTS})
      set(python_path $ENV{PYTHONPATH})
      if(python_path)
        set(python_path ":${python_path}")
      endif(python_path)
      set(python_path "${LIB_STAGE_PATH}/python${PYTHON_VERSION}/site-packages${python_path}")
      if(WIN32)
        set (PY_TESTS_CMD "PYTHONPATH=${python_path}  ${PYTHON_BINARY}")
        # todo fix python unit test running for Windows
        #set (PY_TESTS_CMD "${EXECUTABLE_OUTPUT_PATH}/ost.bat")
        #add_custom_target("${py_test}_run"
        #          CALL "${PY_TESTS_CMD} ${CMAKE_CURRENT_SOURCE_DIR}/${py_test} || echo"
        #          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        #          COMMENT "running checks ${py_test}" VERBATIM)
        #set_target_properties("${py_test}_run" PROPERTIES EXCLUDE_FROM_ALL "1")
        #add_dependencies("${py_test}_run" ost_scripts "_${_ARG_PREFIX}_${_ARG_MODULE}")
        #add_dependencies(check "${py_test}_run")
      else()
        set (PY_TESTS_CMD "PYTHONPATH=${python_path}  ${PYTHON_BINARY}")
        add_custom_target("${py_test}_run"
                  sh -c "${PY_TESTS_CMD} ${CMAKE_CURRENT_SOURCE_DIR}/${py_test}"
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                  COMMENT "running checks ${py_test}" VERBATIM)
        add_custom_target("${py_test}_run_xml"
                  sh -c "${PY_TESTS_CMD} ${CMAKE_CURRENT_SOURCE_DIR}/${py_test} xml || echo"
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                  COMMENT "running checks ${py_test}" VERBATIM)
        add_dependencies("${py_test}_run_xml" ost_scripts "_${_ARG_PREFIX}_${_ARG_MODULE}")
        add_dependencies("${py_test}_run" ost_scripts "_${_ARG_PREFIX}_${_ARG_MODULE}")
        add_dependencies(check "${py_test}_run")
        add_dependencies(check_xml "${py_test}_run_xml")
      endif()
      
    endforeach()
    
    # get_target_property(OUTT ${_test_name} GENERATOR_FILE_NAME)
    # message("${OUTT}")
    
    # get_target_property(OUTT check GENERATOR_FILE_NAME)
    # message("${OUTT}")
endmacro()

#-------------------------------------------------------------------------------
# make sure the previously detected Python interpreter has the given module
#-------------------------------------------------------------------------------
macro(ost_find_python_module MODULE)
  if (NOT PYTHON_MODULE_${MODULE})
    message(STATUS "Searching for python module ${MODULE} for ${PYTHON_BINARY}")
    execute_process(COMMAND ${PYTHON_BINARY} -c "import ${MODULE}"
                    OUTPUT_QUIET ERROR_QUIET
                    RESULT_VARIABLE _IMPORT_ERROR)
    if (_IMPORT_ERROR)
      message(FATAL_ERROR "Could not find python module ${MODULE}. Please install it")
    else()
      message(STATUS "Found python module ${MODULE}")
      set("PYTHON_MODULE_${MODULE}" FOUND CACHE STRING "" FORCE)
    endif()
  endif()
endmacro()

#-------------------------------------------------------------------------------
# make sure the previously detected Python interpreter has the given module
# while allowing alternative module names
#-------------------------------------------------------------------------------
macro(ost_find_python_module_alt MODULES)
  set(_PY_MODS "")
  foreach(py_mod ${MODULES})
    if (NOT PYTHON_MODULE_${py_mod})
      set(_PY_MODS "${_PY_MODS} ${py_mod}")
      message(STATUS "Searching for python module ${py_mod} for ${PYTHON_BINARY}")
      execute_process(COMMAND ${PYTHON_BINARY} -c "import ${py_mod}"
                      OUTPUT_QUIET ERROR_QUIET
                      RESULT_VARIABLE _IMPORT_ERROR)
      if (NOT _IMPORT_ERROR)
        message(STATUS "Found python module ${py_mod}")
        set("PYTHON_MODULE_${py_mod}" FOUND CACHE STRING "" FORCE)
        break()
      endif (NOT _IMPORT_ERROR)
    else ()
      message(STATUS "Already found python module ${py_mod}")
      unset(_IMPORT_ERROR)
      break()
    endif (NOT PYTHON_MODULE_${py_mod})
  endforeach(py_mod ${MODULES})
  if (_IMPORT_ERROR)
    message(FATAL_ERROR "Could not find one of python modules ${_PY_MODS}. Please install one of them.")
  endif ()    
endmacro()

#-------------------------------------------------------------------------------
# this macro tries to detect a very common problem during configuration stage:
# it makes sure that boost python is linked against the same version of python 
# that we are linking against. 
#-------------------------------------------------------------------------------
macro(ost_match_boost_python_version)
  include(CheckCXXSourceRuns)
  # this variable may either be a simple library path or list that contains
  # different libraries for different build-options. For example:
  # optimized;<lib1>;debug;<lib2>
  set(_BOOST_PYTHON_LIBRARY ${BOOST_PYTHON_LIBRARIES})
  list(LENGTH _BOOST_PYTHON_LIBRARY _BP_LENGTH)
  if (_BP_LENGTH GREATER 1)
    list(FIND _BOOST_PYTHON_LIBRARY optimized _OPTIMIZED_INDEX)
    if (_OPTIMIZED_INDEX EQUAL -1)
      message(FATAL_ERROR 
              "Error while trying to get path of boost python library")
    endif()
    math(EXPR _LIB_INDEX 1+${_OPTIMIZED_INDEX})
    list(GET _BOOST_PYTHON_LIBRARY ${_LIB_INDEX} _BP_LIB_PATH)
    set(_BOOST_PYTHON_LIBRARY ${_BP_LIB_PATH})
  endif()
  set(CMAKE_REQUIRED_FLAGS "-I${PYTHON_INCLUDE_PATH} -I${Boost_INCLUDE_DIR} ${PYTHON_LIBRARIES} ${_BOOST_PYTHON_LIBRARY}")
  check_cxx_source_runs(
"#include <boost/python.hpp>

using namespace boost::python;
                           
int main( int argc, char ** argv ) {
  try {
    Py_Initialize();

    object main_module((
      handle<>(borrowed(PyImport_AddModule(\"__main__\")))));

    object main_namespace = main_module.attr(\"__dict__\");

    handle<> ignored(( PyRun_String( \"print 'Hello, World'\",
                                     Py_file_input,
                                     main_namespace.ptr(),
                                     main_namespace.ptr() ) ));
  } catch( error_already_set ) {
    PyErr_Print();
    return 1;
  }   
  return 0;
}" PYTHON_VERSION_FOUND_MATCHES_PYTHON_VERSION_BOOST_WAS_COMPILED_WITH)

  if(NOT PYTHON_VERSION_FOUND_MATCHES_PYTHON_VERSION_BOOST_WAS_COMPILED_WITH)
     message(FATAL_ERROR "Python versions mismatch!\n"
             "Make sure the python version for boost match the one you "
             "specified during configuration\n")
  endif()
endmacro()


#-------------------------------------------------------------------------------
# this macro sets up the stage directories
#-------------------------------------------------------------------------------
macro(setup_stage)
  set(STAGE_DIR "${CMAKE_BINARY_DIR}/stage")
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${STAGE_DIR}/bin  )
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${STAGE_DIR}/bin  )
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${STAGE_DIR}/bin  )
  set(HEADER_STAGE_PATH ${STAGE_DIR}/include )
  set(SHARED_DATA_PATH ${STAGE_DIR}/share/openstructure  )

  if (WIN32)
    set(ENV{PATH} ${STAGE_DIR}/bin\;$ENV{PATH})
  endif (WIN32)
  if (UNIX AND NOT APPLE)
    check_architecture()
  endif()
  set (ARCH ${CMAKE_NATIVE_ARCH})
  if ("${ARCH}" MATCHES "64" AND NOT _UBUNTU_LAYOUT)
    set(LIB_DIR lib64  )
    set(LIB_STAGE_PATH "${STAGE_DIR}/lib64"  )
  else()
    set(LIB_DIR lib  )
    set(LIB_STAGE_PATH "${STAGE_DIR}/lib"  )
  endif()
  if (_UBUNTU_LAYOUT)
    set(LIBEXEC_PATH ${LIB_DIR}/openstructure/libexec  )
    set(LIBEXEC_STAGE_PATH ${LIB_STAGE_PATH}/openstructure/libexec  )
  else()
    set(LIBEXEC_PATH libexec/openstructure  )
    set(LIBEXEC_STAGE_PATH ${STAGE_DIR}/libexec/openstructure  )
  endif()
 
  include_directories("${HEADER_STAGE_PATH}")
  link_directories(${LIB_STAGE_PATH})

endmacro()

#-------------------------------------------------------------------------------
# get compiler version
#-------------------------------------------------------------------------------
function(get_compiler_version _OUTPUT_VERSION)
  exec_program(${CMAKE_CXX_COMPILER}
               ARGS ${CMAKE_CXX_COMPILER_ARG1} -dumpfullversion -dumpversion
               OUTPUT_VARIABLE _COMPILER_VERSION
  )
  string(REGEX REPLACE "([0-9])\\.([0-9])(\\.[0-9])?" "\\1\\2"
    _COMPILER_VERSION ${_COMPILER_VERSION})

  set(${_OUTPUT_VERSION} ${_COMPILER_VERSION} PARENT_SCOPE)
endfunction()



macro(setup_compiler_flags)
  if (WIN32)
     # add_definitions(-DBOOST_TEST_INCLUDED)

     add_definitions(-D_USE_MATH_DEFINES -D_CRT_SECURE_NO_DEPRECATE
                     -D_SCL_SECURE_NO_DEPRECATE -DNOMINMAX)
     #add_definitions(-Zc:wchar_t-)   #
    # add_definitions(-MDd -vmg -EHsc -GR)
    #GR:Uses the __fastcall calling convention (x86 only).
    #-EHsc to specify the synchronous exception handling mode/
    #-vmg Uses full generality for pointers to members.
    #add_definitions(-DBOOST_ZLIB_BINARY=zdll)
    #add_definitions(-NODEFAULTLIB:LIBCMTD.lib)
  endif()


  if (CMAKE_COMPILER_IS_GNUCXX)
    get_compiler_version(_GCC_VERSION)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall" )
    if (_GCC_VERSION MATCHES "44")
      # gcc 4.4. is very strict about aliasing rules. the shared_count
      # implementation that is used boost's shared_ptr violates these rules. To
      # silence the warnings and prevent miscompiles, enable
      #  -fno-strict-aliasing
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-strict-aliasing" )
    endif()
    #message(STATUS "GCC VERSION " ${_GCC_VERSION})
    if (_GCC_VERSION LESS "60")
      # for older compilers we need to enable C++11 for Qt5
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    endif()
    if (_GCC_VERSION LESS "45")
      set(OST_GCC_LESS_45 true)
    else()
      set(OST_GCC_LESS_45 false)
    endif()
  endif()
endmacro()
set(_BOOST_MIN_VERSION 1.31)

macro(setup_boost)
  set (Boost_NO_BOOST_CMAKE TRUE)
  # starting with CMake 3.11 we could use the following instead of the foreach
  # find_package(Boost ${_BOOST_MIN_VERSION} COMPONENTS
  #              python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR} REQUIRED)
  # set(BOOST_PYTHON_LIBRARIES ${Boost_LIBRARIES})
  # see https://cmake.org/cmake/help/v3.11/module/FindBoost.html
  foreach(_python_lib_name python${PYTHON_VERSION_MAJOR}${PYTHON_VERSION_MINOR}
                           python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}
                           python${PYTHON_VERSION_MAJOR}
                           python)
    find_package(Boost ${_BOOST_MIN_VERSION} COMPONENTS ${_python_lib_name} QUIET)
    if(Boost_FOUND)
      message(STATUS "Found Boost package: " ${_python_lib_name})
      set(BOOST_PYTHON_LIBRARIES ${Boost_LIBRARIES})
      break()
    else()
      message(STATUS "Boost package not found: " ${_python_lib_name}
                     ". Trying alternative names!")
    endif()
  endforeach(_python_lib_name)
  if(NOT BOOST_PYTHON_LIBRARIES)
    message(FATAL_ERROR "Failed to find any Boost Python library!")
  endif()
  set(Boost_LIBRARIES)
  find_package(Boost ${_BOOST_MIN_VERSION}
               COMPONENTS unit_test_framework REQUIRED)
  set(BOOST_UNIT_TEST_LIBRARIES ${Boost_LIBRARIES})
  set(Boost_LIBRARIES)
  if (ENABLE_STATIC)
    set(Boost_USE_STATIC_LIBS ON)
  endif()
  find_package(Boost ${_BOOST_MIN_VERSION}
               COMPONENTS filesystem system REQUIRED)
  set(BOOST_LIBRARIES ${Boost_LIBRARIES})
  set(Boost_LIBRARIES)
  find_package(Boost ${_BOOST_MIN_VERSION} COMPONENTS iostreams REQUIRED)
  set(BOOST_IOSTREAM_LIBRARIES ${Boost_LIBRARIES})
  set(Boost_LIBRARIES)
  find_package(Boost ${_BOOST_MIN_VERSION} COMPONENTS program_options REQUIRED)
  set(BOOST_PROGRAM_OPTIONS ${Boost_LIBRARIES})
  set(Boost_LIBRARIES)
  find_package(Boost ${_BOOST_MIN_VERSION} COMPONENTS regex REQUIRED)
  set(BOOST_REGEX_LIBRARIES ${Boost_LIBRARIES})
  set(Boost_LIBRARIES)
  find_package(Boost ${_BOOST_MIN_VERSION} COMPONENTS thread REQUIRED)
  set(BOOST_THREAD ${Boost_LIBRARIES})
  set(Boost_LIBRARIES)
endmacro()


#-------------------------------------------------------------------------------
# Synopsis:
#   ost_action_init()
#
# Description:
#   Initialise cached variables
#-------------------------------------------------------------------------------
macro(ost_action_init)
  set(OST_ACTION_NAMES "" CACHE INTERNAL "" FORCE)
endmacro(ost_action_init)

#-------------------------------------------------------------------------------
# Synopsis:
#   ost_action(ACTION TARGET)
#
# Description:
#   Add a script to actions.
#   ACTION script to be added (needs to have permissions to be executed)
#   TARGET make target to add the action to
#-------------------------------------------------------------------------------
macro(ost_action ACTION TARGET)
  copy_if_different("${CMAKE_CURRENT_SOURCE_DIR}"
                    "${STAGE_DIR}/${LIBEXEC_PATH}" 
                    "${ACTION}" "TARGETS" ${TARGET})
  install(FILES "${ACTION}" DESTINATION "${LIBEXEC_PATH}"
          PERMISSIONS WORLD_EXECUTE GROUP_EXECUTE OWNER_EXECUTE 
                      WORLD_READ GROUP_READ OWNER_READ)
  # storing tool names for bash completion
  string(REGEX REPLACE "^ost-" "" stripped_action ${ACTION})
  if(DEFINED OST_ACTION_NAMES)
    if(${OST_ACTION_NAMES} MATCHES "${stripped_action}")
      set(_ACTION_NAMES "${OST_ACTION_NAMES}")
    else()
      if("${OST_ACTION_NAMES}" STREQUAL "")
        set(_ACTION_NAMES "${stripped_action}")
      else()
        set(_ACTION_NAMES "${OST_ACTION_NAMES} ${stripped_action}")
      endif()
    endif()
  else()
    set(_ACTION_NAMES "${stripped_action}")
  endif()
  set(OST_ACTION_NAMES "${_ACTION_NAMES}" CACHE INTERNAL "" FORCE)
endmacro(ost_action)
