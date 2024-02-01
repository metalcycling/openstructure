#-------------------------------------------------------------------------------
#  Author:  Andreas Schenk
#
#  This file contains functions to determine, copy and verify the dependencies
#  of the application bundle.
#
#  Part of this code was taken from BundleUtilities.cmake, distributed
#  under the OSI-approved BSD License (see below), and adapted.
#-------------------------------------------------------------------------------

#  CMake - Cross Platform Makefile Generator
#  Copyright 2000-2011 Kitware, Inc., Insight Software Consortium
#  All rights reserved.
#  
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions
#  are met:
#  
#  * Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#  
#  * Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#  
#  * Neither the names of Kitware, Inc., the Insight Software Consortium,
#    nor the names of their contributors may be used to endorse or promote
#    products derived from this software without specific prior written
#    permission.
#  
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
#  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
#  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
#  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
#  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
#  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
#  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
#  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#=============================================================================
# function copy dependencies
#=============================================================================
function(copy_dependencies binaries)
  #--------------------------------------------------------------------------
  message("copy dependencies: resolving binary paths")
  #--------------------------------------------------------------------------
  message(STATUS "  binaries='${binaries}'")
  set(resolved_binaries)
  foreach(binary ${binaries})
    get_filename_component(resolved_binary ${binary} REALPATH)
    if(NOT ${resolved_binary} MATCHES ${CMAKE_INSTALL_PREFIX})
      message(FATAL_ERROR "Binary path: '${resolved_binary}' points to outside bundle.")
    endif(NOT ${resolved_binary} MATCHES ${CMAKE_INSTALL_PREFIX})
    list(APPEND resolved_binaries ${resolved_binary})
  endforeach(binary)
  list(REMOVE_DUPLICATES resolved_binaries)

  #--------------------------------------------------------------------------
  message("copy dependencies: determining dependencies")
  #--------------------------------------------------------------------------
  collect_dependencies("${resolved_binaries}" keys)

  #--------------------------------------------------------------------------
  message("copy dependencies: copying files")
  #--------------------------------------------------------------------------
  list(LENGTH keys n)
  math(EXPR n ${n}*2)

  set(i 0)
  foreach(key ${keys})
    math(EXPR i ${i}+1)
    if(${${key}_COPYFLAG})
      message(STATUS "${i}/${n}: copying '${${key}_RESOLVED_ITEM}' to '${${key}_RESOLVED_EMBEDDED_ITEM}'")
      set(item "${${key}_ITEM}")
      copy_resolved_item("${${key}_RESOLVED_ITEM}" "${${key}_RESOLVED_EMBEDDED_ITEM}")
    else(${${key}_COPYFLAG})
      message(STATUS "${i}/${n}: *NOT* copying '${${key}_RESOLVED_ITEM}'")
    endif(${${key}_COPYFLAG})
  endforeach(key)

  #--------------------------------------------------------------------------
  message("copy dependencies: changing install names")
  #--------------------------------------------------------------------------
  foreach(key ${keys})
    math(EXPR i ${i}+1)
    if(APPLE)
      if(NOT ${${key}_SYSTEMFLAG})
        message(STATUS "${i}/${n}: changing install names for '${${key}_RESOLVED_EMBEDDED_ITEM}'")
        change_install_names_for_item(keys "${key}")
      else(NOT ${${key}_SYSTEMFLAG})
        message(STATUS "${i}/${n}: skipping install name change for system library: '${${key}_RESOLVED_EMBEDDED_ITEM}'")
      endif(NOT ${${key}_SYSTEMFLAG})
    else(APPLE)
      message(STATUS "${i}/${n}: changing install names not required on this platform '${${key}_RESOLVED_EMBEDDED_ITEM}'")
    endif(APPLE)
  endforeach(key)


  #--------------------------------------------------------------------------
  message("copy dependencies: verifying binaries")
  #--------------------------------------------------------------------------
  verify(keys)

  #--------------------------------------------------------------------------
  message("copy dependencies: cleaning up")
  #--------------------------------------------------------------------------
  clear_keys(keys)

  #--------------------------------------------------------------------------
  message("copy dependencies: done")
  #--------------------------------------------------------------------------
endfunction(copy_dependencies)

#=============================================================================
# function collect_dependencies
#=============================================================================
function(collect_dependencies binaries keys_var)
  set(${keys_var} PARENT_SCOPE)

  # First add all binaries to the keys. 
  # This is done before dependency resolution because one of the binaries could
  # potentially be a dependency. 
  foreach(exe ${binaries})
    get_filename_component(exepath "${exe}" PATH)
    set_keys_for_item(${keys_var} "${exe}" "${exe}" 0 0)
  endforeach(exe)

  # Get dependencies and add them to the keys.
  foreach(exe ${binaries})
    get_filename_component(exepath "${exe}" PATH)
    add_dependencies_for_item(${keys_var}  "${exe}" "${exepath}")
  endforeach(exe)

  # Propagate values to caller's scope:
  set(${keys_var} ${${keys_var}} PARENT_SCOPE)
  foreach(key ${${keys_var}})
    set(${key}_ITEM "${${key}_ITEM}" PARENT_SCOPE)
    set(${key}_RESOLVED_ITEM "${${key}_RESOLVED_ITEM}" PARENT_SCOPE)
    set(${key}_RESOLVED_EMBEDDED_ITEM "${${key}_RESOLVED_EMBEDDED_ITEM}" PARENT_SCOPE)
    set(${key}_COPYFLAG "${${key}_COPYFLAG}" PARENT_SCOPE)
    set(${key}_SYSTEMFLAG "${${key}_SYSTEMFLAG}" PARENT_SCOPE)
  endforeach(key)

endfunction(collect_dependencies)




#=============================================================================
# function set_keys_for_item
#=============================================================================
function(set_keys_for_item keys_var context resolved_item copyflag systemflag)
  get_item_key("${resolved_item}" key)
  list(LENGTH ${keys_var} length_before)
  append_unique(${keys_var} "${key}")
  list(LENGTH ${keys_var} length_after)

  if(NOT length_before EQUAL length_after)
    if(${copyflag} AND (NOT ${systemflag}) )
      resolve_embedded_item("${resolved_item}" resolved_embedded_item)
      set(${key}_RESOLVED_EMBEDDED_ITEM "${resolved_embedded_item}" PARENT_SCOPE)
      if(${resolved_embedded_item} STREQUAL ${resolved_item})
        set(${key}_COPYFLAG 0 PARENT_SCOPE)
      else(${resolved_embedded_item} STREQUAL ${resolved_item})
        set(${key}_COPYFLAG 1 PARENT_SCOPE)
      endif(${resolved_embedded_item} STREQUAL ${resolved_item})
      set(${key}_SYSTEMFLAG 0 PARENT_SCOPE)
    else(${copyflag} AND (NOT ${systemflag}) )
      set(${key}_COPYFLAG ${copyflag} PARENT_SCOPE)
      set(${key}_SYSTEMFLAG ${systemflag} PARENT_SCOPE)
      set(${key}_RESOLVED_EMBEDDED_ITEM "${resolved_item}" PARENT_SCOPE)
    endif(${copyflag} AND (NOT ${systemflag}) )
    set(${keys_var} ${${keys_var}} PARENT_SCOPE)
    get_filename_component(item_name "${resolved_item}" NAME)
    set(${key}_ITEM "${item_name}" PARENT_SCOPE)
    set(${key}_RESOLVED_ITEM "${resolved_item}" PARENT_SCOPE)

  else(NOT length_before EQUAL length_after)
    message(FATAL_ERROR "Item key '${key}' already in the list. This should not happen.")
  endif(NOT length_before EQUAL length_after)
endfunction(set_keys_for_item)

#=============================================================================
# function get_item_key
#=============================================================================
function(get_item_key item key_var)
  get_filename_component(item_name "${item}" REALPATH)
  if(WIN32)
    string(TOLOWER "${item_name}" item_name)
  endif()
  string(REGEX REPLACE "\\." "_" item_name "${item_name}")
  string(REGEX REPLACE "/" "_"   item_name "${item_name}")
  string(REGEX REPLACE "\\\\" "_" ${key_var} "${item_name}")
  set(${key_var} ${${key_var}} PARENT_SCOPE)
endfunction(get_item_key)


#=============================================================================
# function clear_keys
#=============================================================================
function(clear_keys keys_var)
  foreach(key ${${keys_var}})
    set(${key}_ITEM PARENT_SCOPE)
    set(${key}_RESOLVED_ITEM PARENT_SCOPE)
    set(${key}_RESOLVED_EMBEDDED_ITEM PARENT_SCOPE)
    set(${key}_COPYFLAG PARENT_SCOPE)
    set(${key}_SYSTEMFLAG PARENT_SCOPE)
  endforeach(key)
  set(${keys_var} PARENT_SCOPE)
endfunction(clear_keys)



#=============================================================================
# function resolve_item
#=============================================================================
function(resolve_item context item exepath resolved_item_var)
  set(resolved_item "${item}")

  # @rpath references are relative to the paths built into the binaries with -rpath
  # rpath resolution comes first, as the rpath might contain a @loader_path or @executable_path reference
  if(resolved_item MATCHES "@rpath")
    get_rpath("${context}" rpath)
    string(REPLACE "@rpath/" "${rpath}/" resolved_item "${resolved_item}")
  endif(resolved_item MATCHES "@rpath")

  # @executable_path references are assumed relative to exepath
  if(resolved_item MATCHES "@executable_path")
    string(REPLACE "@executable_path" "${exepath}" resolved_item "${resolved_item}")
  endif(resolved_item MATCHES "@executable_path")

  # @loader_path references are assumed relative to the
  # PATH of the given "context" (presumably another library)
  if(resolved_item MATCHES "@loader_path")
    get_filename_component(contextpath "${context}" PATH)
    string(REPLACE "@loader_path" "${contextpath}" resolved_item "${resolved_item}")
  endif(resolved_item MATCHES "@loader_path")
  
  if(NOT WIN32)
    get_filename_component(resolved_item "${resolved_item}" REALPATH)
  endif(NOT WIN32)

  if(IS_ABSOLUTE "${resolved_item}" AND EXISTS "${resolved_item}")
    set(${resolved_item_var} "${resolved_item}" PARENT_SCOPE)
    return()
  endif(IS_ABSOLUTE "${resolved_item}" AND EXISTS "${resolved_item}")

  set(ri "ri-NOTFOUND")
  if(UNIX AND NOT WIN32)
    find_file(ri "${item}" ${exepath}  NO_DEFAULT_PATH)
    find_file(ri "${item}" ${exepath}  /usr/lib)
    if(APPLE)
      find_file(ri "${item}" PATHS ENV DYLD_LIBRARY_PATH  NO_DEFAULT_PATH)
    else(APPLE)
      find_file(ri "${item}" PATHS  ENV LD_LIBRARY_PATH NO_DEFAULT_PATH)
    endif(APPLE)
  else(UNIX AND NOT WIN32)
    find_program(ri "${item}" PATHS "${exepath}" NO_DEFAULT_PATH)
    find_program(ri "${item}" PATHS "${CMAKE_INSTALL_PREFIX}/bin" NO_DEFAULT_PATH)
    find_program(ri "${item}" PATHS "${exepath}")
  endif(UNIX AND NOT WIN32)
  if(ri)
    set(resolved_item "${ri}")
    set(ri "ri-NOTFOUND")
  endif(ri)

  if(IS_ABSOLUTE "${resolved_item}" AND EXISTS "${resolved_item}")
    set(${resolved_item_var} "${resolved_item}" PARENT_SCOPE)
  else(IS_ABSOLUTE "${resolved_item}" AND EXISTS "${resolved_item}")
    if(APPLE)
      # see if item is available in System frameworks
      string(REPLACE "/Library/Frameworks" "/System/Library/Frameworks" system_resolved_item ${resolved_item})
      if(IS_ABSOLUTE "${system_resolved_item}" AND EXISTS "${system_resolved_item}") 
        set(${resolved_item_var} "${system_resolved_item}" PARENT_SCOPE)
      else(IS_ABSOLUTE "${system_resolved_item}" AND EXISTS "${system_resolved_item}")
         message(FATAL_ERROR "cannot resolve item '${item}' in context:'${context}'. Partially resolved path is: '${system_resolved_item}'") 
      endif(IS_ABSOLUTE "${system_resolved_item}" AND EXISTS "${system_resolved_item}") 
    else(APPLE)
      message(FATAL_ERROR "cannot resolve item '${item}' in context:'${context}'. Partially resolved path is: '${resolved_item}'")
    endif(APPLE)
  endif(IS_ABSOLUTE "${resolved_item}" AND EXISTS "${resolved_item}")
endfunction(resolve_item)


#=============================================================================
# function resolve_embedded_item
#=============================================================================
function(resolve_embedded_item item resolved_item_var)
  if(${item} MATCHES ${CMAKE_INSTALL_PREFIX})
    # item already embedded
    set(${resolved_item_var} ${item} PARENT_SCOPE)
  else(${item} MATCHES ${CMAKE_INSTALL_PREFIX})
    # only embed libraries, therefore put into lib dir
    get_filename_component(item_name "${item}" NAME)
    if(WIN32 AND NOT UNIX)
      set(${resolved_item_var} "${CMAKE_INSTALL_PREFIX}/bin/${item_name}" PARENT_SCOPE)
    else(WIN32 AND NOT UNIX)
      set(${resolved_item_var} "${CMAKE_INSTALL_PREFIX}/${LIB_DIR}/${item_name}" PARENT_SCOPE)
    endif(WIN32 AND NOT UNIX)
  endif(${item} MATCHES ${CMAKE_INSTALL_PREFIX})
endfunction(resolve_embedded_item)





#=============================================================================
# function change_install_names_for_item
#=============================================================================
function(change_install_names_for_item keys_var key)

  # Ensure the item is within the install directory
  if(${${key}_RESOLVED_EMBEDDED_ITEM} MATCHES ${CMAKE_INSTALL_PREFIX})
    set(prereqs "")
    get_filename_component(exepath ${${key}_RESOLVED_EMBEDDED_ITEM} PATH)
    string(REPLACE ${CMAKE_INSTALL_PREFIX} "" rel_path ${exepath})
    string(REGEX REPLACE "/[^/]*" "/.." inv_rel_path ${rel_path})
    get_dependencies_for_item("${${key}_RESOLVED_EMBEDDED_ITEM}" prereqs)
    set(changes "")
    foreach(pr ${prereqs})
      # Each referenced item's key is "rkey" in the loop:
      get_filename_component(exepath ${${key}_RESOLVED_EMBEDDED_ITEM} PATH)
      resolve_item( ${${key}_RESOLVED_EMBEDDED_ITEM} ${pr} exepath resolved_pr)
      # If dependency already resolves to library within install directory we don't have to do anything.
      if(NOT (${resolved_pr} MATCHES ${CMAKE_INSTALL_PREFIX}) AND EXISTS "${resolved_pr}" )
        get_item_key("${resolved_pr}" rkey)
        list(FIND ${keys_var} ${rkey} r_found)
        if(r_found EQUAL -1)
          message("Error while trying to get dependencies for ${${key}_RESOLVED_EMBEDDED_ITEM}.")
          MESSAGE(FATAL_ERROR "Key: '${rkey}' not found in dependency list. Cannot change install name.")
        endif(r_found EQUAL -1)
        if(NOT ${${rkey}_SYSTEMFLAG})
          get_filename_component(ref_exepath ${${rkey}_RESOLVED_EMBEDDED_ITEM} PATH)
          string(REPLACE ${CMAKE_INSTALL_PREFIX} "" ref_rel_path ${ref_exepath})
          set(diff_path @loader_path${inv_rel_path}${ref_rel_path}/)
          set(changes ${changes} "-change" "${pr}" "${diff_path}${${rkey}_ITEM}")
        endif(NOT ${${rkey}_SYSTEMFLAG})
      endif(NOT (${resolved_pr} MATCHES ${CMAKE_INSTALL_PREFIX}) AND EXISTS "${resolved_pr}" )
    endforeach(pr)
     
    execute_process(COMMAND chmod u+w "${${key}_RESOLVED_EMBEDDED_ITEM}")
    
    # Change this item's id and all of its references in one call
    # to install_name_tool:
    execute_process(COMMAND install_name_tool ${changes} -id "${${key}_ITEM}" "${${key}_RESOLVED_EMBEDDED_ITEM}")
  else( ${${key}_RESOLVED_EMBEDDED_ITEM} MATCHES ${CMAKE_INSTALL_PREFIX})
    MESSAGE(FATAL_ERROR "Trying to change install name for: ${${key}_RESOLVED_EMBEDDED_ITEM}, which is outside of install directory")
  endif(${${key}_RESOLVED_EMBEDDED_ITEM} MATCHES ${CMAKE_INSTALL_PREFIX})

endfunction(change_install_names_for_item)



#=============================================================================
# function verify
#=============================================================================
function(verify keys_var)
  foreach(key ${${keys_var}})
    if(NOT ${${key}_SYSTEMFLAG})
      get_dependencies_for_item("${${key}_RESOLVED_EMBEDDED_ITEM}" dependencies)
      get_filename_component(exepath ${${key}_RESOLVED_EMBEDDED_ITEM} PATH)
      foreach(dep ${dependencies})
        resolve_item( ${${key}_RESOLVED_EMBEDDED_ITEM} ${dep} exepath resolved_dep)
        is_system_lib(${resolved_dep} system_flag)
        if(NOT ${system_flag})
          if(NOT ${resolved_dep} MATCHES ${CMAKE_INSTALL_PREFIX})
            MESSAGE("Warning: item:'${${key}_RESOLVED_EMBEDDED_ITEM}' contains external dependency:'${resolved_dep}'")
          endif(NOT ${resolved_dep} MATCHES ${CMAKE_INSTALL_PREFIX})
        endif(NOT ${system_flag})
      endforeach(dep)
    endif(NOT ${${key}_SYSTEMFLAG})
  endforeach(key)
endfunction(verify)


#=============================================================================
# function add_dependencies_for_item
#=============================================================================
function(add_dependencies_for_item keys_var  item exepath)
  if(NOT IS_ABSOLUTE "${item}" AND EXISTS "${item}")
    message(FATAL_ERROR "Item '${item}' does not exist.")
  endif(NOT IS_ABSOLUTE "${item}" AND EXISTS "${item}")

  get_dependencies_for_item(${item} dependencies)
  foreach(dep ${dependencies})
    resolve_item("${item}" "${dep}" "${exepath}" resolved_dep)
    get_item_key("${resolved_dep}" dep_key)
    list(FIND ${keys_var} ${dep_key} dep_found)
    if(dep_found EQUAL -1)
      message(STATUS "Resolved '${dep}' to: '${resolved_dep}'")
      is_system_lib(${resolved_dep} system_flag)
      if(${system_flag})
        MESSAGE(STATUS "Skipping system library:'${resolved_dep}'")
        set_keys_for_item(${keys_var} ${item} ${resolved_dep} 0 1)
      else(${system_flag})
        set_keys_for_item(${keys_var} ${item} ${resolved_dep} 1   0)
        add_dependencies_for_item(${keys_var} ${resolved_dep} ${exepath})
      endif(${system_flag})
    endif(dep_found EQUAL -1)
  endforeach(dep ${dependencies})
  # Propagate values to caller's scope:
  set(${keys_var} ${${keys_var}} PARENT_SCOPE)
  foreach(key ${${keys_var}})
    set(${key}_ITEM "${${key}_ITEM}" PARENT_SCOPE)
    set(${key}_RESOLVED_ITEM "${${key}_RESOLVED_ITEM}" PARENT_SCOPE)
    set(${key}_RESOLVED_EMBEDDED_ITEM "${${key}_RESOLVED_EMBEDDED_ITEM}" PARENT_SCOPE)
    set(${key}_COPYFLAG "${${key}_COPYFLAG}" PARENT_SCOPE)
    set(${key}_SYSTEMFLAG "${${key}_SYSTEMFLAG}" PARENT_SCOPE)
  endforeach(key)
endfunction(add_dependencies_for_item)



#=============================================================================
# function append_unique
#=============================================================================
function(append_unique list_var value)
  set(contains 0)
  foreach(item ${${list_var}})
    if("${item}" STREQUAL "${value}")
      set(contains 1)
      break()
    endif("${item}" STREQUAL "${value}")
  endforeach(item)

  if(NOT contains)
    set(${list_var} ${${list_var}} "${value}" PARENT_SCOPE)
  endif(NOT contains)
endfunction(append_unique)



#=============================================================================
# function read_config
#=============================================================================
function(read_config ost_config_var)
  file(READ "${CMAKE_INSTALL_PREFIX}/${LIBEXEC_PATH}/ost_config" ost_config)
  set(${ost_config_var} "${ost_config}" PARENT_SCOPE)
endfunction(read_config)

#=============================================================================
# function write_config
#=============================================================================
function(write_config ost_config)
  file(WRITE "${CMAKE_INSTALL_PREFIX}/${LIBEXEC_PATH}/ost_config" "${ost_config}")
endfunction(write_config)

#=============================================================================
# function write_bundle_binary
#=============================================================================
function(write_bundle_binary binary_name real_binary)
  file(WRITE "${CMAKE_INSTALL_PREFIX}/../MacOS/${binary_name}" "#!/bin/sh\nexec \"\${0%/*}/${real_binary}\"\n")
  execute_process(COMMAND chmod ugo+x "${CMAKE_INSTALL_PREFIX}/../MacOS/${binary_name}")
endfunction(write_bundle_binary)

#=============================================================================
# function add_linked_bundle
#=============================================================================
function(add_linked_bundle parent_bundle_name bundle_name binary_name bundle_identifier bundle_icon)
  get_filename_component(dmg_main_path ${CMAKE_INSTALL_PREFIX}/../../../../ PATH)
  file(MAKE_DIRECTORY "${dmg_main_path}/${bundle_name}.app/Contents/MacOS/")
  file(MAKE_DIRECTORY "${dmg_main_path}/${bundle_name}.app/Contents/Resources/")
  file(WRITE "${dmg_main_path}/${bundle_name}.app/Contents/MacOS/${binary_name}" "#!/bin/sh\nexec \"\${0%/*}/../../../${parent_bundle_name}.app/Contents/Resources/bin/${binary_name}\"\n")
  get_filename_component(icon_name ${bundle_icon} NAME)
  file(COPY ${bundle_icon} DESTINATION "${dmg_main_path}/${bundle_name}.app/Contents/Resources")
  create_plist(${binary_name} ${bundle_name} ${bundle_identifier} ${icon_name} "${dmg_main_path}/${bundle_name}.app/Contents/Info.plist")
  execute_process(COMMAND chmod ugo+x "${dmg_main_path}/${bundle_name}.app/Contents/MacOS/${binary_name}")
endfunction(add_linked_bundle)

#=============================================================================
# function create_plist
#=============================================================================
function(create_plist BUNDLE_EXECUTABLE BUNDLE_NAME BUNDLE_IDENTIFIER BUNDLE_ICON FILENAME)
  set( info_plist "<?xml version=\"1.0\" encoding=\"UTF-8\"?>
<!DOCTYPE plist PUBLIC \"-//Apple Computer//DTD PLIST 1.0//EN\" 
    \"http://www.apple.com/DTDs/PropertyList-1.0.dtd\">
<plist version=\"1.0\">
<dict>
  <key>CFBundleDevelopmentRegion</key>
  <string>English</string>
  <key>CFBundleExecutable</key>
  <string>\@BUNDLE_EXECUTABLE\@</string>
  <key>CFBundleName</key>
  <string>\@BUNDLE_NAME\@</string>  
  <key>CFBundleDisplayName</key>
  <string>\@BUNDLE_NAME\@</string>  
  <key>CFBundleHelpBookFolder</key>
  <string>\@BUNDLE_NAME\@ Manual</string>
  <key>CFBundleHelpBookName</key>
  <string>OpenStructure Help</string>
  <key>CFBundleIconFile</key>
  <string>\@BUNDLE_ICON\@</string>
  <key>CFBundleIdentifier</key>
  <string>\@BUNDLE_IDENTIFIER\@</string>
  <key>CFBundleInfoDictionaryVersion</key>
  <string>6.0</string>
  <key>CFBundlePackageType</key>
  <string>APPL</string>
  <key>CFBundleSignature</key>
  <string>????</string>
  <key>CFBundleVersion</key>
  <string></string>
  <key>NSMainNibFile</key>
  <string>MainMenu</string>
  <key>NSPrincipalClass</key>
  <string>NSApplication</string>
</dict>
</plist>")
  string(CONFIGURE "${info_plist}" info_plist_replaced)
  file(WRITE  ${FILENAME} "${info_plist_replaced}")
endfunction(create_plist)

#=============================================================================
# function echo_all_cmake_variable_values
#=============================================================================

function(echo_all_cmake_variable_values)
  message("")
  get_cmake_property(vs VARIABLES)
  foreach(v ${vs})
    message("${v}='${${v}}'")
  endforeach(v)
  message("")
endfunction(echo_all_cmake_variable_values)


#=============================================================================
#=============================================================================
# OSX specific functions
#=============================================================================
#=============================================================================

if(APPLE)

#=============================================================================
# function get_dependencies_for_item (OSX)
#=============================================================================
function(get_dependencies_for_item  item list_var)

  if(NOT IS_ABSOLUTE "${item}" AND EXISTS "${item}")
    message(FATAL_ERROR "Cannot find item: '${item}' to get dependencies.")
  endif(NOT IS_ABSOLUTE "${item}" AND EXISTS "${item}")

  find_program(otool_cmd "otool" PATHS "/usr/local/bin" "/usr/bin")
  if(NOT otool_cmd)
    message(FATAL_ERROR "Could not find otool - cannot analyze dependencies.")
  endif(NOT otool_cmd)

  execute_process( COMMAND ${otool_cmd} -L ${item} OUTPUT_VARIABLE otool_cmd_ov RESULT_VARIABLE retcode )
  if(retcode)
    MESSAGE(FATAL_ERROR "otool stopped with return code: '${retcode}'")
  endif(retcode)
  # Convert to a list of lines:
  string(REGEX REPLACE ";" "\\\\;" candidates "${otool_cmd_ov}")
  set(eol_char "E")
  string(REGEX REPLACE "\n" "${eol_char};" candidates "${candidates}")

  # check for install id and remove it from list, since otool -L can include a reference to itself
  set(install_id)
  execute_process( COMMAND ${otool_cmd} -D ${item} OUTPUT_VARIABLE install_id_ov  RESULT_VARIABLE retcode)
  if(retcode)
    MESSAGE(FATAL_ERROR "otool stopped with return code: '${retcode}'")
  endif(retcode)
  # second line is install name
  string(REGEX REPLACE ".*:\n" "" install_id "${install_id_ov}")
  if(install_id)
    string(REGEX MATCH "[^\n ].*[^\n ]" install_id "${install_id}")
  endif(install_id)

  # Analyze each line for file names that match the regular expression:
  set(otool_regex "^\t([^\t]+) \\(compatibility version ([0-9]+.[0-9]+.[0-9]+), current version ([0-9]+.[0-9]+.[0-9]+)\\)${eol_char}$")
  set(dep_list)
  foreach(candidate ${candidates})
    if("${candidate}" MATCHES "${otool_regex}")
      # Extract information from each candidate:
      string(REGEX REPLACE "${otool_regex}" "\\1" dep "${candidate}")
      if(NOT "${dep}" STREQUAL "${install_id}")
        list(APPEND dep_list ${dep})
      endif(NOT "${dep}" STREQUAL "${install_id}")
    endif("${candidate}" MATCHES "${otool_regex}")
  endforeach(candidate)
  set(${list_var} ${dep_list} PARENT_SCOPE)
endfunction(get_dependencies_for_item)

#=============================================================================
# function get_rpath (OSX)
#=============================================================================
function(get_rpath target rpath_var)
  find_program(rpath_exe "otool" PATHS "/usr/local/bin" "/usr/bin")
  if(NOT rpath_exe)
    message(FATAL_ERROR "Could not find otool - cannot resolve rpath.")
  endif(NOT rpath_exe)
 execute_process(
    COMMAND ${rpath_exe} "-l" ${target}
    OUTPUT_VARIABLE rpath_cmd_ov
    RESULT_VARIABLE retcode
  )
  if(retcode)
    MESSAGE(FATAL_ERROR "otool stopped with return code: '${retcode}'")
  endif(retcode)
  set(rpath_cmd_regex ".*LC_RPATH\n *cmdsize [0-9]*\n\ *path ([^\n ]*) \\(offset.*")
  if(rpath_cmd_ov MATCHES ${rpath_cmd_regex})
    string(REGEX REPLACE ${rpath_cmd_regex} "\\1"  rpath ${rpath_cmd_ov})
    set(${rpath_var} ${rpath} PARENT_SCOPE)
  endif(rpath_cmd_ov MATCHES ${rpath_cmd_regex})
endfunction(get_rpath)

#=============================================================================
# function copy_qt (OSX)
#=============================================================================
function(copy_qt library_dir plugin_dir plugins)
  find_file(QT_MENU_NIB qt_menu.nib PATHS ${library_dir}/QtGui.framework/Versions/4 ${library_dir} PATH_SUFFIXES Resources)
  if(NOT QT_MENU_NIB)
    message(FATAL "qt_menu.nib not found")
  endif(NOT QT_MENU_NIB)
  file(COPY "${QT_MENU_NIB}"
       DESTINATION "${CMAKE_INSTALL_PREFIX}/${LIBEXEC_PATH}/"
       FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
  file(COPY ${plugin_dir}
       DESTINATION ${CMAKE_INSTALL_PREFIX}/${LIBEXEC_PATH}) 
  file(GLOB_RECURSE QT_PLUGINS "${CMAKE_INSTALL_PREFIX}/${LIBEXEC_PATH}/*.dylib")
  set(${plugins} ${QT_PLUGINS} PARENT_SCOPE)
  file(WRITE "${CMAKE_INSTALL_PREFIX}/${LIBEXEC_PATH}/qt.conf" "[Paths]\nPlugins=../plugins\n")
endfunction(copy_qt)

#=============================================================================
# function is_system_lib (OSX)
#=============================================================================
function(is_system_lib item system_var)
  if(item MATCHES "^(/System/Library/|/usr/lib/)")
    set(${system_var}  1 PARENT_SCOPE)
  else()
    set(${system_var}  0 PARENT_SCOPE)
  endif()
endfunction(is_system_lib)

#=============================================================================
# function copy_resolved_item (OSX)
#=============================================================================
function(copy_resolved_item resolved_item resolved_embedded_item)
  if("${resolved_item}" STREQUAL "${resolved_embedded_item}")
    message(STATUS "warning: resolved_item == resolved_embedded_item - not copying...")
  else()
    message(STATUS "copying COMMAND ${CMAKE_COMMAND} -E copy ${resolved_item} ${resolved_embedded_item}")
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy "${resolved_item}" "${resolved_embedded_item}")
  endif()

endfunction(copy_resolved_item)

#=============================================================================
# function copy_python (OSX)
#=============================================================================
function(copy_python include_path version new_binary_path)


  get_filename_component(real_python_include_path ${include_path} REALPATH)
  get_filename_component(python_root_dir ${real_python_include_path}/../.. REALPATH) 
  file(COPY ${python_root_dir}/Resources/Python.app/Contents/MacOS/Python DESTINATION ${CMAKE_INSTALL_PREFIX})
  file(COPY ${python_root_dir}/${LIB_DIR}/python${version}/ DESTINATION ${CMAKE_INSTALL_PREFIX}/${LIB_DIR}/python${version})
  file(GLOB  py_config_files "${include_path}/pyconfig*.h")
  file(COPY ${py_config_files} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/python${version})
  file(GLOB_RECURSE python_so_files "${CMAKE_INSTALL_PREFIX}/${LIB_DIR}/python${version}/*.so")
  foreach(so_file ${python_so_files})
    remove_rpath( "${so_file}")
  endforeach(so_file)
  read_config(ost_config)
  file(COPY ${python_root_dir}/Resources/Python.app/Contents/MacOS/Python DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
  set(${new_binary_path} "${CMAKE_INSTALL_PREFIX}/bin/Python" PARENT_SCOPE)
  string(REGEX REPLACE "pyexec=\"[^\n\$]*\"" "pyexec=\"\$DNG_BINDIR/Python\"\nexport PYTHONHOME=\$DNG_ROOT" ost_config "${ost_config}")
  write_config("${ost_config}")
endfunction(copy_python)

#=============================================================================
# function remove_rpath (OSX)
#=============================================================================
function(remove_rpath file)
  get_rpath(${file} rpath)
  if(${rpath})
    execute_process(COMMAND stat -f %Mp%Lp "${file}" OUTPUT_VARIABLE permissions)
    execute_process(COMMAND chmod ugo+wx "${file}")
    string(SUBSTRING ${permissions} 0 3 permissions)
    file(RPATH_REMOVE FILE "${file}")
    execute_process(COMMAND chmod ${permissions} "${file}")
  endif(${rpath})
endfunction(remove_rpath)

elseif(WIN32 AND NOT UNIX)

#=============================================================================
#=============================================================================
# Windows specific functions
#=============================================================================
#=============================================================================


#=============================================================================
# function get_dependencies_for_item (Windows)
#=============================================================================
function(get_dependencies_for_item  item list_var)

  if(NOT IS_ABSOLUTE "${item}" AND EXISTS "${item}")
    message(FATAL_ERROR "Cannot find item: '${item}' to get dependencies.")
  endif(NOT IS_ABSOLUTE "${item}" AND EXISTS "${item}")

  find_program(dumpbin_cmd "dumpbin")
  if(NOT dumpbin_cmd)
    message(FATAL_ERROR "Could not find dumpbin - cannot analyze dependencies.")
    return()
  endif(NOT dumpbin_cmd)

  set(ENV{VS_UNICODE_OUTPUT} "") # Block extra output from inside VS IDE.
  execute_process( COMMAND ${dumpbin_cmd}  /dependents ${item} OUTPUT_VARIABLE dumpbin_cmd_ov  ERROR_VARIABLE dumpbin_cmd_ev RESULT_VARIABLE retcode)
  if(retcode)
    if(${item} MATCHES "System32")
      # 32 bit cmake might run within WOW64, which automatically mangles the path to System32
      # dumpbin, which is x86_64 in this case will run with the unmangled path and therefore has to look for
      # the 32 bit libraries in SysWOW64 instead of System32
      string(REPLACE "System32" "SysWOW64" item ${item})
      execute_process( COMMAND ${dumpbin_cmd}  /dependents ${item} OUTPUT_VARIABLE dumpbin_cmd_ov  ERROR_VARIABLE dumpbin_cmd_ev RESULT_VARIABLE retcode)
    endif(${item} MATCHES "System32")
  endif(retcode)
  if(retcode)
    MESSAGE( "dumpbin output: '${dumpbin_cmd_ov}'")
    MESSAGE( "dumpbin error output: '${dumpbin_cmd_ev}'")
    MESSAGE(FATAL_ERROR "dumpbin stopped with return code: '${retcode}'")
  endif(retcode)
  # Convert to a list of lines:
  string(REGEX REPLACE ";" "\\\\;" candidates "${dumpbin_cmd_ov}")
  set(eol_char "E")
  string(REGEX REPLACE "\n" "${eol_char};" candidates "${candidates}")

  # Analyze each line for file names that match the regular expression:
  set(dumpbin_regex "^    ([^ ].*[Dd][Ll][Ll])${eol_char}$")
  set(dep_list)
  foreach(candidate ${candidates})
    if("${candidate}" MATCHES "${dumpbin_regex}")
      # Extract information from each candidate:
      string(REGEX REPLACE "${dumpbin_regex}" "\\1" dep "${candidate}")
      list(APPEND dep_list ${dep})
    endif("${candidate}" MATCHES "${dumpbin_regex}")
  endforeach(candidate)
  set(${list_var} ${dep_list} PARENT_SCOPE)
endfunction(get_dependencies_for_item)

#=============================================================================
# function get_rpath (Windows)
#=============================================================================
function(get_rpath target rpath_var)
  message(FATAL_ERROR "get_rpath not implemented for Windows.")
endfunction(get_rpath)

#=============================================================================
# function is_system_lib (Windows)
#=============================================================================
function(is_system_lib item system_var)
  string(TOLOWER "${item}" lower)
  string(TOLOWER "$ENV{SystemRoot}" sysroot)
  string(REGEX REPLACE "\\\\" "/" sysroot "${sysroot}")
  string(TOLOWER "$ENV{windir}" windir)
  string(REGEX REPLACE "\\\\" "/" windir "${windir}")
  if(lower MATCHES "^(${sysroot}/sys(tem|wow)|${windir}/sys(tem|wow)|(.*/)*msvc[^/]+dll)" AND NOT lower MATCHES "python.*dll")
    set(${system_var}  1 PARENT_SCOPE)
  else()
    set(${system_var}  0 PARENT_SCOPE)
  endif()
endfunction(is_system_lib)

#=============================================================================
# function copy_resolved_item (Windows)
#=============================================================================
function(copy_resolved_item resolved_item resolved_embedded_item)
  # ignore case on Windows
  string(TOLOWER "${resolved_item}" resolved_item_compare)
  string(TOLOWER "${resolved_embedded_item}" resolved_embedded_item_compare)

  if("${resolved_item_compare}" STREQUAL "${resolved_embedded_item_compare}")
    message(STATUS "warning: resolved_item == resolved_embedded_item - not copying...")
  else()
    message(STATUS "copying COMMAND ${CMAKE_COMMAND} -E copy ${resolved_item} ${resolved_embedded_item}")
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy "${resolved_item}" "${resolved_embedded_item}")
  endif()

endfunction(copy_resolved_item)

#=============================================================================
# function copy_python (Windows)
#=============================================================================
function(copy_python include_path version new_binary_path)
  get_filename_component(real_python_include_path ${include_path} REALPATH)
  get_filename_component(python_root_dir ${real_python_include_path}/.. REALPATH) 
  file(GLOB  lib_files "${python_root_dir}/Lib/*")
  file(COPY ${lib_files} DESTINATION ${CMAKE_INSTALL_PREFIX}/${LIB_DIR}/python${version})
  file(GLOB  dll_files "${python_root_dir}/DLLs/*")
  file(COPY ${dll_files} DESTINATION ${CMAKE_INSTALL_PREFIX}/${LIB_DIR}/python${version})
  file(GLOB  py_config_files "${include_path}/pyconfig*.h")
  file(COPY ${py_config_files} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/python${version})
  file(COPY ${python_root_dir}/python.exe DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
  set(${new_binary_path} "${CMAKE_INSTALL_PREFIX}/bin/python.exe" PARENT_SCOPE)
  
  # copy msvc90 CRT for python
  get_filename_component(msvc90_dir "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VisualStudio\\9.0;InstallDir]" ABSOLUTE)
  get_filename_component(msvc90_express_dir "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\VCExpress\\9.0;InstallDir]" ABSOLUTE)
  find_path(msvc_redist_dir NAMES x86/Microsoft.VC90.CRT/Microsoft.VC90.CRT.manifest
            PATHS
            "${msvc90_dir}/../../VC/redist"
            "${msvc90_express_dir}/../../VC/redist"
            "C:/Program Files/Microsoft Visual Studio 9.0/VC/redist"
            "C:/Program Files (x86)/Microsoft Visual Studio 9.0/VC/redist"
    )
  set(msvc90_crt_dir "${msvc_redist_dir}/x86/Microsoft.VC90.CRT")
  file(COPY "${msvc90_crt_dir}/Microsoft.VC90.CRT.manifest"
            "${msvc90_crt_dir}/msvcm90.dll"
            "${msvc90_crt_dir}/msvcp90.dll"
            "${msvc90_crt_dir}/msvcr90.dll"
       DESTINATION "${CMAKE_INSTALL_PREFIX}/bin")

endfunction(copy_python)

#=============================================================================
# function copy_qt (Windows)
#=============================================================================
function(copy_qt library_dir plugin_dir plugins)
  file(COPY ${plugin_dir}
       DESTINATION ${CMAKE_INSTALL_PREFIX}/${LIBEXEC_PATH}) 
  file(GLOB_RECURSE QT_PLUGINS "${CMAKE_INSTALL_PREFIX}/${LIBEXEC_PATH}-/*.dll")
  set(${plugins} ${QT_PLUGINS} PARENT_SCOPE)
  file(WRITE "${CMAKE_INSTALL_PREFIX}/${LIBEXEC_PATH}/qt.conf" "[Paths]\nPlugins=../plugins\n")
endfunction(copy_qt)


elseif(UNIX)


#=============================================================================
#=============================================================================
# Linux specific functions
#=============================================================================
#=============================================================================

#=============================================================================
# function get_dependencies_for_item (Linux)
#=============================================================================
function(get_dependencies_for_item item list_var)

  if(NOT IS_ABSOLUTE "${item}" AND EXISTS "${item}")
    message(FATAL_ERROR "Cannot find item: '${item}' to get dependencies.")
  endif(NOT IS_ABSOLUTE "${item}" AND EXISTS "${item}")

  find_program(ldd_cmd "ldd" PATHS "/usr/local/bin" "/usr/bin")
  if(NOT ldd_cmd)
    message(FATAL_ERROR "Could not find ldd - cannot analyze dependencies.")
    return()
  endif(NOT ldd_cmd)
  set(old_ld_library_path $ENV{LD_LIBRARY_PATH})
  set(ENV{LD_LIBRARY_PATH} "${CMAKE_INSTALL_PREFIX}/${LIB_DIR}:$ENV{LD_LIBRARY_PATH}")
  execute_process( COMMAND ${ldd_cmd}  ${item} OUTPUT_VARIABLE ldd_cmd_ov ERROR_VARIABLE ldd_cmd_ev RESULT_VARIABLE retcode)
  set(ENV{LD_LIBRARY_PATH} ${old_ld_library_path})
  if(retcode)
    if(${ldd_cmd_ov})
      message(${ldd_cmd_ov})
    endif(${ldd_cmd_ov})
    if(${ldd_cmd_ev})
      message(${ldd_cmd_ev})
    endif(${ldd_cmd_ev})
    message("command:'${ldd_cmd}'  '${item}'")
    message(FATAL_ERROR "ldd stopped with return code: '${retcode}'")
  endif(retcode)
  # Convert to a list of lines:
  string(REGEX REPLACE ";" "\\\\;" candidates "${ldd_cmd_ov}")
  set(eol_char "E")
  string(REGEX REPLACE "\n" "${eol_char};" candidates "${candidates}")

  # Analyze each line for file names that match the regular expression:
  set(ldd_regex "^[\t ]*[^\t ]+ => ([^\t\(]+) .*${eol_char}$")
  set(ldd_regex_error "not found${eol_char}$")
  set(ldd_regex_fallback "^[\t ]*([^\t ]+) => ([^\t ]+).*${eol_char}$")
  set(dep_list)
  foreach(candidate ${candidates})
    if("${candidate}" MATCHES "${ldd_regex}")
      # Extract information from each candidate:
      if("${candidate}" MATCHES "${ldd_regex_error}")
        string(REGEX REPLACE "${ldd_regex_fallback}" "\\1" dep "${candidate}")
      else("${candidate}" MATCHES "${ldd_regex_error}")
        string(REGEX REPLACE "${ldd_regex}" "\\1" dep "${candidate}")
      endif("${candidate}" MATCHES "${ldd_regex_error}")
      list(APPEND dep_list ${dep})
    endif("${candidate}" MATCHES "${ldd_regex}")
  endforeach(candidate)
  set(${list_var} ${dep_list} PARENT_SCOPE)
endfunction(get_dependencies_for_item)

#=============================================================================
# function get_rpath (Linux)
#=============================================================================
function(get_rpath target rpath_var)
  find_program(rpath_exe "readelf" PATHS "/usr/local/bin" "/usr/bin")
  if(NOT rpath_exe)
    message(FATAL_ERROR "Could not find readelf - cannot resolve rpath.")
  endif(NOT rpath_exe)
 execute_process(
    COMMAND ${rpath_exe} "-d" ${target}
    OUTPUT_VARIABLE rpath_cmd_ov
    RESULT_VARIABLE retcode
  )
  if(retcode)
    MESSAGE(FATAL_ERROR "readelf stopped with return code: '${retcode}'")
  endif(retcode)
  
  set(rpath_cmd_regex "Library rpath: \[([^\n ]*):\]")
  if(rpath_cmd_ov MATCHES ${rpath_cmd_regex})
    string(REGEX REPLACE ${rpath_cmd_regex} "\\1"  rpath ${rpath_cmd_ov})
    set(${rpath_var} ${rpath} PARENT_SCOPE)
  endif(rpath_cmd_ov MATCHES ${rpath_cmd_regex})
endfunction(get_rpath)

#=============================================================================
# function copy_qt (Linux)
#=============================================================================
function(copy_qt library_dir plugin_dir plugins)
  file(COPY ${plugin_dir}
       DESTINATION ${CMAKE_INSTALL_PREFIX}/${LIBEXEC_PATH}) 
  file(GLOB_RECURSE QT_PLUGINS "${CMAKE_INSTALL_PREFIX}/${LIBEXEC_PATH}/*.so")
  foreach(plugin ${QT_PLUGINS})
    remove_rpath("${plugin}")
  endforeach(plugin)
  set(${plugins} ${QT_PLUGINS} PARENT_SCOPE)
  file(WRITE "${CMAKE_INSTALL_PREFIX}/${LIBEXEC_PATH}/qt.conf" "[Paths]\nPlugins=../plugins\n")
endfunction(copy_qt)

#=============================================================================
# function is_system_lib (Linux)
#=============================================================================
function(is_system_lib item system_var)
  if(item MATCHES "/cmov/|ld-linux[\\.\\-]|libc[\\.\\-]|libdl[\\.\\-]|libdrm[\\.\\-]|libelf[\\.\\-]|libexpat[\\.\\-]|libfontconfig[\\.\\-]|libfreetype[\\.\\-]|libg[\\.\\-]|libgcc_s[\\.\\-]|libGL[\\.\\-]|libglib[\\.\\-]|libgthread[\\.\\-]|lib(ice|ICE)[\\.\\-]|libnvidia[\\.\\-]|libpthread[\\.\\-]|libse(pol|linux)[\\.\\-]|libSM[\\.\\-]|libstdc\\+\\+[\\.\\-]|libX[a-zA-Z0-9]*[\\.\\-]|libxcb[\\.\\-]|pylibc[\\.\\-]")
    set(${system_var} 1 PARENT_SCOPE)
  else()
    set(${system_var}  0 PARENT_SCOPE)
  endif()
endfunction(is_system_lib)

#=============================================================================
# function copy_resolved_item (Linux)
#=============================================================================
function(copy_resolved_item resolved_item resolved_embedded_item)
  get_filename_component(real_resolved_item "${resolved_item}" REALPATH)
  get_filename_component(real_resolved_embedded_item "${resolved_embedded_item}" REALPATH)

  if("${real_resolved_item}" STREQUAL "${real_resolved_embedded_item}")
    message(STATUS "warning: resolved_item == resolved_embedded_item - not copying...")
  else()
    message(STATUS "copying ${resolved_item} to ${resolved_embedded_item} including symlinks")
    # determine all the symlinks pointing to the current item and copy them into the package
    # along with the item
    get_filename_component(external_path  ${real_resolved_item} PATH)
    get_filename_component(embedded_path  ${real_resolved_embedded_item} PATH)
    get_filename_component(external_name  ${real_resolved_item} NAME)
    string(REGEX REPLACE "([a-zA-Z]+).*" "\\1" name_base ${external_name})
    file(GLOB candidates "${external_path}/${name_base}*" )
    set(external_files)
    foreach(candidate ${candidates})
      get_filename_component(real_candidate "${candidate}" REALPATH)
      if(real_candidate STREQUAL ${real_resolved_item})
        list(APPEND external_files ${candidate})
      endif(real_candidate STREQUAL ${real_resolved_item})
    endforeach(candidate)
    file(COPY ${external_files} DESTINATION ${embedded_path})
    remove_rpath("${real_resolved_embedded_item}")
  endif()

endfunction(copy_resolved_item)

#=============================================================================
# function copy_python (Linux)
#=============================================================================
function(copy_python include_path version new_binary_path)
  if (_UBUNTU_LAYOUT)
    if ("${ARCH}" MATCHES "64")
      set(UBUNTU_ARCH_STRING "x86_64")
    else()  
      set(UBUNTU_ARCH_STRING "i386")
    endif()
  endif()
  get_filename_component(real_python_include_path ${include_path} REALPATH)
  get_filename_component(python_root_dir ${real_python_include_path}/../.. REALPATH) 
  if (_UBUNTU_LAYOUT)
    file (GLOB_RECURSE python_lib_files RELATIVE ${python_root_dir}/${LIB_DIR}/python${version}/ FOLLOW_SYMLINKS ${python_root_dir}/${LIB_DIR}/python${version}/*)
    foreach(python_lib_file ${python_lib_files})
      if(NOT ${python_lib_file} MATCHES "pyc$") 
        get_filename_component(real_python_lib_file ${python_root_dir}/${LIB_DIR}/python${version}/${python_lib_file} REALPATH)
        execute_process(COMMAND cmake -E copy ${real_python_lib_file} ${CMAKE_INSTALL_PREFIX}/${LIB_DIR}/python${version}/${python_lib_file})
      endif()
    endforeach()
  else()
    file(COPY ${python_root_dir}/${LIB_DIR}/python${version}/ DESTINATION ${CMAKE_INSTALL_PREFIX}/${LIB_DIR}/python${version} PATTERN "*.pyc" EXCLUDE)
    if("${LIB_DIR}" MATCHES "lib64")
      # for 64 bit systems we have to copy the python dir within lib in addition as norarch site packages might get installed there
      file(COPY ${python_root_dir}/lib/python${version}/ DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/python${version} PATTERN "*.pyc" EXCLUDE)
    endif("${LIB_DIR}" MATCHES "lib64")
  endif()
  file(GLOB  py_config_files "${include_path}/pyconfig*.h")
  file(COPY ${py_config_files} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/python${version})
  file(GLOB_RECURSE python_so_files "${CMAKE_INSTALL_PREFIX}/${LIB_DIR}/python${version}/*.so")
  read_config(ost_config)
  if (_UBUNTU_LAYOUT)
    file(GLOB python_libs "${python_root_dir}/${LIB_DIR}/python${version}/config-${UBUNTU_ARCH_STRING}-linux-gnu/libpython${version}.so.*" FOLLOW_SYMLINKS)
  else()
    file(GLOB python_libs "${python_root_dir}/${LIB_DIR}/libpython${version}.so*")
  endif()
  file(COPY ${python_libs} DESTINATION ${CMAKE_INSTALL_PREFIX}/${LIB_DIR})
  foreach( additional_directory $ENV{PYTHONPATH})
    file(COPY ${additional_directory}/ DESTINATION ${CMAKE_INSTALL_PREFIX}/${LIB_DIR}/python${version}/site-packages PATTERN "*.pyc" EXCLUDE)
  endforeach()
  file(COPY ${python_root_dir}/bin/python${version} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
  set(${new_binary_path} "${CMAKE_INSTALL_PREFIX}/bin/python${version}" PARENT_SCOPE)
  string(REGEX REPLACE "pyexec=\"[^\n\$]*\"" "pyexec=\"\$DNG_BINDIR/python${version}\"\nexport PYTHONHOME=\$DNG_ROOT" ost_config "${ost_config}")
  foreach(so_file ${python_so_files})
    file(RPATH_REMOVE FILE "${so_file}")
  endforeach(so_file)
  write_config("${ost_config}")
endfunction(copy_python)


#=============================================================================
# function remove_rpath (Linux)
#=============================================================================
function(remove_rpath file)
  execute_process(COMMAND stat --format %a  "${file}" OUTPUT_VARIABLE permissions)
  execute_process(COMMAND chmod ugo+x "${file}")
  string(SUBSTRING ${permissions} 0 3 permissions)
  file(RPATH_REMOVE FILE "${file}")
  execute_process(COMMAND chmod ${permissions} "${file}")
endfunction(remove_rpath)



else(APPLE)
  MESSAGE("Unknown platform, deployment will not work.")
endif(APPLE)
