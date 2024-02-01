# - helper module to find OSX frameworks

IF(NOT CMAKE_FIND_FRAMEWORKS_INCLUDED)
  SET(CMAKE_FIND_FRAMEWORKS_INCLUDED 1)
  MACRO(CMAKE_FIND_FRAMEWORKS fwk)
    SET(${fwk}_FRAMEWORKS)
    IF(APPLE)
      FOREACH(dir
          /opt/local/Library/Frameworks
          ~/Library/Frameworks
          /Library/Frameworks
          /System/Library/Frameworks
          /Network/Library/Frameworks)
        IF(EXISTS "${dir}/${fwk}.framework")
          SET(${fwk}_FRAMEWORKS "${dir}/${fwk}.framework")
          SET(${fwk}_INCLUDE_PATH "${dir}/${fwk}.framework/Headers")
          SET(${fwk}_LIBRARY_PATH "${dir}") 
          BREAK()
        ENDIF()
      ENDFOREACH(dir)
    ENDIF(APPLE)
  ENDMACRO(CMAKE_FIND_FRAMEWORKS)
ENDIF(NOT CMAKE_FIND_FRAMEWORKS_INCLUDED)
