## This file should be placed in the root directory of your project.
## Then modify the CMakeLists.txt file in the root directory of your
## project to incorporate the testing dashboard.
## # The following are required to uses Dart and the Cdash dashboard
##   ENABLE_TESTING()
##   INCLUDE(CTest)
set(CTEST_PROJECT_NAME "OpenStructure")
set(CTEST_NIGHTLY_START_TIME "00:00:00 MEZ")

set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "kb-pdwg1059-06.pz.unibas.ch")
set(CTEST_DROP_LOCATION "/cdash/submit.php?project=OpenStructure")
set(CTEST_DROP_SITE_CDASH TRUE)

if ("${CMAKE_NATIVE_ARCH}" MATCHES "64")
  set(CMAKE_SYSTEM_NAME "${CMAKE_SYSTEM_NAME}64bit")
endif()
if (USE_DOUBLE_PRECISION)
  set(CMAKE_SYSTEM_NAME "${CMAKE_SYSTEM_NAME}-double")
endif()
if (NOT ENABLE_GFX)
  set(CMAKE_SYSTEM_NAME "${CMAKE_SYSTEM_NAME}-no_gfx")
endif()
if (NOT ENABLE_UI)
  set(CMAKE_SYSTEM_NAME "${CMAKE_SYSTEM_NAME}-no_ui")
endif()

set(BUILDNAME "${CMAKE_SYSTEM_NAME}")
