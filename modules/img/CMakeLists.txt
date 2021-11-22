add_subdirectory(base)
add_subdirectory(alg)
if (IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/ex)
  add_subdirectory(ex)
endif()

