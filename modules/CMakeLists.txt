# recursively process all subdirectories
file(GLOB ALL_FILES *)
foreach(full_name ${ALL_FILES})
  get_filename_component(fname ${full_name} NAME)
  if(NOT ${fname} MATCHES ".*\\.svn" 
     AND NOT ${fname} MATCHES ".*CMakeFiles" 
     AND NOT ${fname} MATCHES ".*main"
     AND NOT ${fname} MATCHES "doc")
      if(IS_DIRECTORY ${full_name})
        add_subdirectory(${fname})
      endif()
  endif()
endforeach()

