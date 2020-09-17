#-------------------------------------------------------------------------------
# Author: Marco Biasini
#-------------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.10.2 FATAL_ERROR)

# Set CMake policies
# Behaviour of target_link_libraries, always link by full path. CMP0060 can
# disable this to stay compatible with old projects. We want the new behaviour
# since before we needed to provide the full path ourselves, that is: we
# manually forced CMake to the new scheme.
cmake_policy(SET CMP0060 NEW)

project(OpenStructure CXX C)
set (CMAKE_EXPORT_COMPILE_COMMANDS 1)
set (OST_VERSION_MAJOR 2)
set (OST_VERSION_MINOR 1)
set (OST_VERSION_PATCH 0)
set (OST_VERSION_STRING ${OST_VERSION_MAJOR}.${OST_VERSION_MINOR}.${OST_VERSION_PATCH} )
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake_support)
include(OST)

option(USE_SHADER "whether to compile with shader support"
       OFF)
option(USE_RPATH "embed rpath upon make install"
       OFF)
option(COMPILE_TMTOOLS "whether to compile the tmalign and tmscore programs"
       OFF)
option(PROFILE "whether to compile with profiling support"
       OFF)       
option(ENABLE_GUI "whether the graphical user interface should be enabled"
       ON)
option(ENABLE_GFX "whether graphics support should be enabled"
       ON)
option(OPENGLPREFERENCE_LEGACY "whether to use the OpenGL Legacy GL implementation"
       OFF)
option(ENABLE_IMG "whether the image processing module should be compiled"
       ON)
option(ENABLE_INFO "whether openstructure should be compiled with support for the info library"
       ON)
option(USE_NUMPY "whether numpy support is added"
      OFF)
option(USE_DOUBLE_PRECISION "whether to compile in double precision" 
       OFF)
option(ENABLE_SPNAV "whether 3DConnexion devices should be supported"
      OFF)
option(STATIC_PROPERTY_WORKAROUND "workaround for static property bug with some boost/boost_python combinations" OFF)
option(COMPILE_TESTS "whether unit tests should be compiled by default" OFF)
option(ENABLE_STATIC "whether static libraries should be compiled" OFF)
option(UBUNTU_LAYOUT "whether Debian/Ubuntu's lib and libexec directory layout should be used" OFF)
option(HIDDEN_VISIBILITY "on gcc, use -fvisibility=hidden" OFF)
option(ENABLE_MM "whether openmm support is added" OFF)

if (CXX)
  set(CMAKE_CXX_COMPILER ${CXX})
endif()

if (CC)
  set(CMAKE_C_COMPILER ${CC})
endif()

if (FORTRAN_COMPILER)
  set(CMAKE_FORTRAN_COMPILER ${FORTRAN_COMPILER})
endif()



if (PREFIX)
  set(CMAKE_INSTALL_PREFIX ${PREFIX})
endif()

if (SYS_ROOT)
  set(CMAKE_PREFIX_PATH ${SYS_ROOT})
endif()

if (OPTIMIZE)
  set(CMAKE_BUILD_TYPE Release)
  set(_OPT ON)
else()
   set(CMAKE_BUILD_TYPE Debug)
    set(_OPT OFF)
endif()

if (USE_DOUBLE_PRECISION)
  set(_DOUBLE_PREC ON)
else()
  set(_DOUBLE_PREC OFF)
  set(_DOUBLE_PREC OFF)
endif()

if (NOT ENABLE_INFO)
  set(ENABLE_GFX OFF)
endif()

if (NOT ENABLE_GFX)
  set(ENABLE_GUI OFF)
endif()

if (ENABLE_INFO)
  set(_INFO ON)
else()
  set(_INFO OFF)
endif()

if (ENABLE_GUI)
  set(_UI ON)
else()
  set(_UI OFF)
endif()

if (ENABLE_GFX)
  set(_OPENGL ON)
else()
  set(_OPENGL OFF)
endif()

if (USE_SHADER)
  set(_SHADER ON)
else()
  set(_SHADER OFF)
endif()
if (USE_NUMPY)
  set(_NUMPY ON)
else()
  set(_NUMPY OFF)
endif()
if (COMPILE_TMTOOLS)
  set(_TM_TOOLS ON)
else()
  set(_TM_TOOLS OFF)
endif()
if (ENABLE_IMG)
  set(_IMG ON)
else()
  set(_IMG OFF)
endif()
if (PROFILE)
  set(_PROFILE ON)
else()
  set(_PROFILE OFF)
endif()
if(ENABLE_MM)
  set(_OPENMM ON)
else()
  set(_OPENMM OFF)
endif()


if (UBUNTU_LAYOUT)
  set(_UBUNTU_LAYOUT ON)
else()
  set(_UBUNTU_LAYOUT OFF)
endif()

if (COMPOUND_LIB)
  set(_COMP_LIB "${COMPOUND_LIB}")
  if (NOT IS_ABSOLUTE "${COMPOUND_LIB}")
    set(COMPOUND_LIB "${CMAKE_CURRENT_BINARY_DIR}/${COMPOUND_LIB}")
  endif()
else()
  set(_COMP_LIB "NONE")
endif()

if (STATIC_PROPERTY_WORKAROUND)
  set(_STATICPROPS ON)
else()
  set(_STATICPROPS OFF)
endif()

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()

setup_stage()
file(MAKE_DIRECTORY ${STAGE_DIR}
     ${EXECUTABLE_OUTPUT_PATH}
     ${HEADER_STAGE_PATH}
     ${LIB_STAGE_PATH}
     ${LIBEXEC_STAGE_PATH})

setup_compiler_flags()

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY)

if (ENABLE_STATIC AND UNIX AND NOT APPLE)
  set (CMAKE_FIND_LIBRARY_SUFFIXES .a)
endif()

# Python needed before Boost
find_package(Python 3.6 REQUIRED)
# Split version string
string(REPLACE "." ";" _python_version_list ${PYTHON_VERSION})
list(GET _python_version_list 0 PYTHON_VERSION_MAJOR)
list(GET _python_version_list 1 PYTHON_VERSION_MINOR)

if(BOOST_VERSION)
  message(STATUS "forcing boost to version ${BOOST_VERSION}")
  set(_BOOST_MIN_VERSION ${BOOST_VERSION})
endif()
setup_boost()
if (WIN32)
   ADD_DEFINITIONS("-DBOOST_ALL_NO_LIB")
endif()

if (ENABLE_INFO)        
  find_package(Qt5Xml REQUIRED)
  set(QT_LIBRARIES ${QT_LIBRARIES} Qt5::Xml)
endif()

if (OPTIMIZE)
  if (CXX_FLAGS)
    set (CMAKE_CXX_FLAGS_RELEASE ${CXX_FLAGS})
  endif()
  if (C_FLAGS)
    set (CMAKE_C_FLAGS_RELEASE ${C_FLAGS})
  endif()
else()
  if (CXX_FLAGS)
    set(CMAKE_CXX_FLAGS_DEBUG ${CXX_FLAGS})
  endif()
  if (C_FLAGS)
    set(CMAKE_C_FLAGS_DEBUG ${C_FLAGS})
  endif()
endif()
if (ENABLE_GFX)
  if (OPENGLPREFERENCE_LEGACY)
    set(OpenGL_GL_PREFERENCE LEGACY)
  else()
    set(OpenGL_GL_PREFERENCE GLVND)
  endif(OPENGLPREFERENCE_LEGACY)
  find_package(OpenGL REQUIRED)
endif()

if (ENABLE_IMG OR ENABLE_GFX)
  find_package(PNG REQUIRED)
else()
  set (PNG_LIBRARY "")
endif()

find_package(Eigen 3.2.0 REQUIRED)
find_package(SQLite3 3.7.13 REQUIRED)

if(USE_NUMPY)
  find_package(Numpy REQUIRED)
endif()

if(ENABLE_MM)
  find_package(OpenMM REQUIRED)
  set(_OPENMM_PLUGINS "${OPEN_MM_PLUGIN_DIR}")
else(ENABLE_MM)
  set(_OPENMM_PLUGINS "NONE")
endif(ENABLE_MM)

if (ENABLE_IMG)
  find_package(FFTW REQUIRED)
  find_package(TIFF REQUIRED)
endif()


if (ENABLE_STATIC)
  set(Boost_LIBRARIES)
  set(Boost_USE_STATIC_LIBS ON)
  find_package(Boost ${_BOOST_MIN_VERSION} 
               COMPONENTS filesystem system iostreams regex thread REQUIRED)
  find_package(ZLIB REQUIRED)
   
  if (UNIX AND NOT APPLE)
    set(STATIC_LIBRARIES ${Boost_LIBRARIES} ${FFTW_LIBRARIES} ${TIFF_LIBRARIES} 
        ${PNG_LIBRARY} ${ZLIB_LIBRARY})
  else ()
    set(STATIC_LIBRARIES ${Boost_LIBRARIES} ${FFTW_LIBRARIES} ${TIFF_LIBRARIES} 
        ${PNG_LIBRARY} ${ZLIB_LIBRARY})
  endif()
endif()

if (ENABLE_GUI)
  find_package(Qt5Widgets REQUIRED)
  find_package(Qt5OpenGL REQUIRED)
  set(QT_LIBRARIES ${QT_LIBRARIES} Qt5::Gui Qt5::Widgets Qt5::Core)
  ost_find_python_module(PyQt5)
  ost_find_python_module_alt("PyQt5.sip;sip")
endif()

if (ENABLE_SPNAV)
  #find_package(SpNav REQUIRED)
  set(_SPNAV ON)
else()
  set(_SPNAV OFF)
endif()

if (CMAKE_COMPILER_IS_GNUCXX)
  # do not write back into cache, otherwise the compile command line gets expanded
  # with multiple -fno-strict-aliasing flags, triggering a complete rebuild whenever
  # cmake is run
  set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fno-strict-aliasing")
  if("${CMAKE_CXX_COMPILER_VERSION}" VERSION_GREATER "4.6")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-attributes")
  endif("${CMAKE_CXX_COMPILER_VERSION}" VERSION_GREATER "4.6")
  if(HIDDEN_VISIBILITY) 
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden")
    set (_HIDDEN_VIS ON)
  else(HIDDEN_VISIBILITY)
    set (_HIDDEN_VIS OFF)
  endif(HIDDEN_VISIBILITY)
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-warning-option")
endif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
if(CMAKE_C_COMPILER_ID MATCHES "Clang")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-unknown-warning-option")
endif(CMAKE_C_COMPILER_ID MATCHES "Clang")

# basic environment
include_directories(${Boost_INCLUDE_DIRS} 
                    ${FFTW_INCLUDE_DIRS} 
                    ${EIGEN3_INCLUDE_DIRS}
                    ${TIFF_INCLUDE_DIR}
                    ${PNG_INCLUDE_DIRS}
                    ${OPEN_MM_INCLUDE_DIRS}
                    ${SQLITE3_INCLUDE_DIRS}
                    )
if(USE_NUMPY)
include_directories(${PYTHON_NUMPY_INCLUDE_DIR})
endif()

if (UNIX)
  SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
  SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
  if(USE_RPATH)
    set(_USE_RPATH ON)
    SET(CMAKE_INSTALL_RPATH "$ORIGIN/../${LIB_DIR}")
    SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  else()
    set(_USE_RPATH OFF)
   SET(CMAKE_INSTALL_RPATH "")
   SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH FALSE)
  endif() 
endif()

if (ENABLE_STATIC)
  SET(CMAKE_SKIP_BUILD_RPATH TRUE)
  SET(LINK_SEARCH_END_STATIC TRUE) 
endif()

add_subdirectory(modules)
add_subdirectory(scripts)
add_subdirectory(tools)
add_subdirectory(actions)
# deployment has to come last, to ensure that all install commands are run before deployment
# magic is done
add_subdirectory(deployment)

set(FILES_TO_BE_REMOVED ${PROJECT_BINARY_DIR}/stage tests)
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
                             "${FILES_TO_BE_REMOVED}")


#ost_match_boost_python_version(${PYTHON_LIBRARIES})

if (CMAKE_COMPILER_IS_GNUCXX)
   set(HIDDEN_VIS_MSG 
       "\n   Hidden object visibility  (-DHIDDEN_VISIBILITY) : ${_HIDDEN_VIS}")
endif()
message(STATUS 
        "OpenStructure will be built with the following options:\n"
        "   Install Prefix                        (-DPREFIX) : ${CMAKE_INSTALL_PREFIX}\n"
        "   RPath in install                   (-DUSE_RPATH) : ${_USE_RPATH}\n"
        "   Info support                     (-DENABLE_INFO) : ${_INFO}\n"
        "   Graphical interface               (-DENABLE_GUI) : ${_UI}\n"
        "   OpenGL support                    (-DENABLE_GFX) : ${_OPENGL}\n"
        "   Image Processing support          (-DENABLE_IMG) : ${_IMG}\n"
        "   Shader support                    (-DUSE_SHADER) : ${_SHADER}\n"
        "   Numpy support                      (-DUSE_NUMPY) : ${_NUMPY}\n"
        "   SpaceNav Device support         (-DENABLE_SPNAV) : ${_SPNAV}\n"
        "   OpenMM support                     (-DENABLE_MM) : ${_OPENMM}\n"
        "   OpenMM plugins            (-DOPEN_MM_PLUGIN_DIR) : ${_OPENMM_PLUGINS}\n"
        "   Optimize                            (-DOPTIMIZE) : ${_OPT}\n"
        "   Profiling support                    (-DPROFILE) : ${_PROFILE}\n"
        "   Double Precision        (-DUSE_DOUBLE_PRECISION) : ${_DOUBLE_PREC}\n"
        "   Compound Lib                    (-DCOMPOUND_LIB) : ${_COMP_LIB}\n"
        "   TMAlign and TMScore          (-DCOMPILE_TMTOOLS) : ${_TM_TOOLS}\n"
        "   Static Libraries               (-DENABLE_STATIC) : ${ENABLE_STATIC}\n"
        "   Debian/Ubuntu directory layout (-DUBUNTU_LAYOUT) : ${_UBUNTU_LAYOUT}"
        ${HIDDEN_VIS_MSG})


# LocalWords:  CMake CMP
