#===============================================================================
# CMake script for Thomas Relativistic Electronic Structure Calculation (TRESC)
# System: Linux / Windows

# ---------------release---------------
# cmake -B build/release -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=ifx
# cmake --build build/release
# chmod +x tshell.sh && cp tshell.sh build
# ----------------debug----------------
# cmake -B build/debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=ifx
# cmake --build build/debug
# chmod +x tshell.sh && cp tshell.sh build
#===============================================================================


cmake_minimum_required(VERSION 3.12)

project(TRESC
  LANGUAGES Fortran
  DESCRIPTION "Thomas Relativistic Electronic Structure Calculation"
)


#-------------------------------------------------------------------------------
# General build settings
#-------------------------------------------------------------------------------
set(CMAKE_Fortran_COMPILER_PREPROCESS ON)
enable_language(Fortran)
set(CMAKE_Fortran_LINKER_LANGUAGE 
"Fortran" CACHE STRING "Set linker language to Fortran" FORCE)
set(CMAKE_MAKE_PROGRAM "ninja" CACHE INTERNAL "")
set(
  CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIR} CACHE PATH "" FORCE
)

# Disable C/C++ toolchain autodetection
set(CMAKE_C_COMPILER_WORKS 1)
set(CMAKE_CXX_COMPILER_WORKS 1)

set(CMAKE_Fortran_STANDARD 2018)
set(CMAKE_Fortran_MODULE_DIR ${CMAKE_BINARY_DIR}/modules)
file(MAKE_DIRECTORY ${CMAKE_Fortran_MODULE_DIR})
set(CMAKE_VERBOSE_MAKEFILE ON)

# Ninja generator preferences
set(CMAKE_SUPPRESS_REGENERATION true)
set(CMAKE_SKIP_INSTALL_ALL_DEPENDENCY true)

set(CMAKE_Fortran_OBJECT_OUTPUT_DIR "${CMAKE_BINARY_DIR}")
file(MAKE_DIRECTORY ${CMAKE_Fortran_OBJECT_OUTPUT_DIR})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/build)


#-------------------------------------------------------------------------------
# Source files and executable
#-------------------------------------------------------------------------------
set(TRESC_SOURCES
  src/Fundamentals.f90
  src/Atoms.f90
  src/Functional.f90
  src/GRysroot.f90
  src/Rysdata.f90
  src/Rys.f90
  src/Lebedev.f90
  src/Representation.f90
  src/Hamiltonian.f90
  src/SCF.f90
  src/Tkernel.f90
)
# check exists
foreach(src_file IN LISTS TRESC_SOURCES)
  if(NOT EXISTS ${CMAKE_SOURCE_DIR}/${src_file})
    message(FATAL_ERROR "Missing source file: ${CMAKE_SOURCE_DIR}/${src_file}")
  endif()
endforeach()
add_executable(TRESC ${TRESC_SOURCES})
set_target_properties(TRESC PROPERTIES
  OUTPUT_NAME_DEBUG "tkerneld"
  OUTPUT_NAME_RELEASE "tkernel"
)
# dependencies
set_source_files_properties(
  ${CMAKE_SOURCE_DIR}/src/Fundamentals.f90
  PROPERTIES
  OBJECT_OUTPUT_NAME "Fundamentals"
  Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIR}
)
set_source_files_properties(
  ${CMAKE_SOURCE_DIR}/src/Atoms.f90
  PROPERTIES
  OBJECT_OUTPUT_NAME "Atoms"
  Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIR}
  OBJECT_DEPENDS ${CMAKE_SOURCE_DIR}/src/Fundamentals.f90
)
set_source_files_properties(
  ${CMAKE_SOURCE_DIR}/src/Rysdata.f90
  PROPERTIES
  OBJECT_OUTPUT_NAME "Rysdata"
  Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIR}
  OBJECT_DEPENDS ${CMAKE_SOURCE_DIR}/src/Fundamentals.f90
)
set_source_files_properties(
  ${CMAKE_SOURCE_DIR}/src/GRysroot.f90
  PROPERTIES
  OBJECT_OUTPUT_NAME "GRysroot"
  Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIR}
  OBJECT_DEPENDS ${CMAKE_SOURCE_DIR}/src/Fundamentals.f90
)
set_source_files_properties(
  ${CMAKE_SOURCE_DIR}/src/Rys.f90
  PROPERTIES
  OBJECT_OUTPUT_NAME "Rys"
  Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIR}
  OBJECT_DEPENDS ${CMAKE_SOURCE_DIR}/src/Rysdata.f90
)
set_source_files_properties(
  ${CMAKE_SOURCE_DIR}/src/Lebedev.f90
  PROPERTIES
  OBJECT_OUTPUT_NAME "Lebedev"
  Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIR}
  OBJECT_DEPENDS ${CMAKE_SOURCE_DIR}/src/Atoms.f90
)
set_source_files_properties(
  ${CMAKE_SOURCE_DIR}/src/Functional.f90
  PROPERTIES
  OBJECT_OUTPUT_NAME "Functional"
  Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIR}
  OBJECT_DEPENDS "${CMAKE_SOURCE_DIR}/src/Atoms.f90"
)
set_source_files_properties(
  ${CMAKE_SOURCE_DIR}/src/Representation.f90
  PROPERTIES
  OBJECT_OUTPUT_NAME "Representation"
  Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIR}
  OBJECT_DEPENDS "${CMAKE_SOURCE_DIR}/src/Lebedev.f90;
${CMAKE_SOURCE_DIR}/src/Functional.f90"
)
set_source_files_properties(
  ${CMAKE_SOURCE_DIR}/src/Hamiltonian.f90
  PROPERTIES
  OBJECT_OUTPUT_NAME "Hamiltonian"
  Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIR}
  OBJECT_DEPENDS "${CMAKE_SOURCE_DIR}/src/GRysroot.f90;
${CMAKE_SOURCE_DIR}/src/Rys.f90;${CMAKE_SOURCE_DIR}/src/Atoms.f90"
)
set_source_files_properties(
  ${CMAKE_SOURCE_DIR}/src/SCF.f90
  PROPERTIES
  OBJECT_OUTPUT_NAME "SCF"
  Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIR}
  OBJECT_DEPENDS "${CMAKE_SOURCE_DIR}/src/Hamiltonian.f90;
${CMAKE_SOURCE_DIR}/src/Representation.f90"
)
set_source_files_properties(
  ${CMAKE_SOURCE_DIR}/src/Tkernel.f90
  PROPERTIES
  OBJECT_OUTPUT_NAME "Tkernel"
  Fortran_MODULE_DIRECTORY ${CMAKE_Fortran_MODULE_DIR}
  OBJECT_DEPENDS ${CMAKE_SOURCE_DIR}/src/SCF.f90
)


#-------------------------------------------------------------------------------
# Dependency handling
#-------------------------------------------------------------------------------
# Check envar
if(NOT DEFINED ENV{ONEAPI_ROOT})
  message(FATAL_ERROR "ONEAPI_ROOT environment variable not set")
else()
  set(ONEAPI_ROOT $ENV{ONEAPI_ROOT})
endif()
if(NOT DEFINED ENV{LIBXC_ROOT})
  message(FATAL_ERROR "LIBXC_ROOT environment variable not set")
else()
  set(LIBXC_ROOT $ENV{LIBXC_ROOT})
endif()
# Find MKL
find_library(MKL_BLAS95_LIB
  NAMES 
    mkl_blas95_lp64
    libmkl_blas95_lp64.a
  PATHS 
    "${ONEAPI_ROOT}/mkl/latest/lib"
)
if(NOT MKL_BLAS95_LIB)
    message(FATAL_ERROR "MKL BLAS95 library not found")
else()
    message(STATUS "Found MKL BLAS95 library: ${MKL_BLAS95_LIB}")
endif()

find_library(MKL_LAPACK95_LIB
  NAMES 
    mkl_lapack95_lp64
    libmkl_lapack95_lp64.a
  PATHS 
    "${ONEAPI_ROOT}/mkl/latest/lib"
)
if(NOT MKL_LAPACK95_LIB)
    message(FATAL_ERROR "MKL LAPACK95 library not found")
else()
    message(STATUS "Found MKL LAPACK95 library: ${MKL_LAPACK95_LIB}")
endif()
# Find LibXC
find_library(XC_LIB
  NAMES 
    xc
    libxc.a
  PATHS 
    "${LIBXC_ROOT}/build"
)
if(NOT XC_LIB)
    message(FATAL_ERROR "xc library not found")
else()
    message(STATUS "Found xc library: ${XC_LIB}")
endif()

find_library(XCF03_LIB
  NAMES 
    xcf03
    libxcf03.a
  PATHS 
    "${LIBXC_ROOT}/build"
)
if(NOT XCF03_LIB)
    message(FATAL_ERROR "xcf03 library not found")
else()
    message(STATUS "Found xcf03 library: ${XCF03_LIB}")
endif()


#-------------------------------------------------------------------------------
# Include directories and linking
#-------------------------------------------------------------------------------
target_include_directories(TRESC PRIVATE
  "${LIBXC_ROOT}/build"
  "${CMAKE_SOURCE_DIR}/src"
  "${ONEAPI_ROOT}/mkl/latest/include/mkl/intel64/ilp64"
  "${ONEAPI_ROOT}/mpi/latest/include"
  "${ONEAPI_ROOT}/mpi/latest/include/mpi/ilp64"
)
set(MKL_INCLUDE_DIR
  "${ONEAPI_ROOT}/mkl/latest/include/mkl/intel64/ilp64"
)
target_compile_options(TRESC PRIVATE
  $<$<PLATFORM_ID:Windows>:/I${LIBXC_ROOT}/build>
  $<$<PLATFORM_ID:Windows>:/I${MKL_INCLUDE_DIR}>
  $<$<NOT:$<PLATFORM_ID:Windows>>:-I${LIBXC_ROOT}/build>
  $<$<NOT:$<PLATFORM_ID:Windows>>:-I${MKL_INCLUDE_DIR}>
  # MPI will be searched automatically
)
target_link_libraries(TRESC PRIVATE
  "${MKL_BLAS95_LIB}"
  "${MKL_LAPACK95_LIB}"
  "${XC_LIB}"
  "${XCF03_LIB}"
)



#-------------------------------------------------------------------------------
# Platform-agnostic compiler configuration
#-------------------------------------------------------------------------------
if(WIN32)
  set(CMAKE_Fortran_COMPILER "ifx" CACHE STRING "Intel Fortran" FORCE)
  # general compile options
  string(CONCAT COMPILE_STR
    "/module:${CMAKE_Fortran_MODULE_DIR} "
    "/nologo "
    "/libs:dll "      # dynamic link to Fortran libs
    "/threads "
    "/assume:byterecl "
    "/Qopenmp "
    "/Qmkl:parallel "
    "/QxCORE-AVX2 "   # /QxHost only works on intel processors
    "/Qfma "
  )
  set(CMAKE_Fortran_FLAGS "${COMPILE_STR}")
  # debug compile options
  string(CONCAT DEBUG_STR
    "/stand:f08 "
    "/dbglibs "
    "/Od "
    "/debug:full "
    "/check:all "
    "/warn:all "
    "/traceback "
  )
  set(CMAKE_Fortran_FLAGS_DEBUG "${DEBUG_STR}")
  # release compile options
  string(CONCAT RELEASE_STR
    "/O3 "
    "/Qipo "
    #"/fp:precise "   # I want more precession
    #"/fp:fast=2 "    # I want more efficiency
  )
  set(CMAKE_Fortran_FLAGS_RELEASE "${RELEASE_STR}")
elseif(UNIX)
  set(CMAKE_Fortran_COMPILER "ifx" CACHE STRING "Intel Fortran" FORCE)
  # general compile options
  string(CONCAT COMPILE_STR
    "-module ${CMAKE_Fortran_MODULE_DIR} "
    "-nologo "
    "-qmkl=parallel "
    "-threads "
    "-assume byterecl "
    "-qopenmp "
    "-xCORE-AVX2 "    # -xHost only works on intel processors
    "-fma "
  )
  set(CMAKE_Fortran_FLAGS "${COMPILE_STR}")
  # debug compile options
  string(CONCAT DEBUG_STR
    "-stand f08 "
    "-O0 "
    "-g "
    "-check all "
    "-warn all "
    "-traceback "
  )
  set(CMAKE_Fortran_FLAGS_DEBUG "${DEBUG_STR}")
  # release compile options
  string(CONCAT RELEASE_STR
    "-O3 "
    "-ipo "
    #"-fp-model precise "   # I want more precession
    #"-fp-model fast=2 "    # I want more efficiency
  )
  set(CMAKE_Fortran_FLAGS_RELEASE "${RELEASE_STR}")
endif()


#-------------------------------------------------------------------------------
# Compile and Link options
#-------------------------------------------------------------------------------
set(CMAKE_Fortran_OUTPUT_EXTENSION ".obj")
if(WIN32)
  string(CONCAT COMPILE_OBJECT_STR
    "<CMAKE_Fortran_COMPILER> <FLAGS> /c <SOURCE> "
    "/o ${CMAKE_Fortran_OBJECT_OUTPUT_DIR}/<OBJECT>"
  )
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    string(CONCAT LINK_EXECUTABLE_STR
      "<CMAKE_Fortran_COMPILER> "
      "/nologo "
      "-fuse-ld=lld "            # link by lld-link
      "/debug "
      "/MDd "
      "/libs:dll "
      "/Qopenmp "
      "/Qmkl:parallel "
      "<OBJECTS> "
      "/o <TARGET> "
      "/link <LINK_LIBRARIES>"   # manually pass to lld-link
    )
  elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    string(CONCAT LINK_EXECUTABLE_STR
      "<CMAKE_Fortran_COMPILER> "
      "/nologo "
      "/MD "
      "/libs:dll "
      "/Qipo "                   # with -fuse-ld=lld
      "/Qopenmp "
      "/Qmkl:parallel "
      "<OBJECTS> "
      "/o <TARGET> "
      "/link <LINK_LIBRARIES>"   # manually pass to lld-link
    )
  endif()
elseif(UNIX)
  string(CONCAT COMPILE_OBJECT_STR
    "<CMAKE_Fortran_COMPILER> <FLAGS> -c <SOURCE> "
    "-o ${CMAKE_Fortran_OBJECT_OUTPUT_DIR}/<OBJECT>"
  )
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    string(CONCAT LINK_EXECUTABLE_STR
      "<CMAKE_Fortran_COMPILER> "
      "-nologo "
      "-fuse-ld=lld "            # link by ld.lld
      "-g "
      "-qopenmp "
      "-qmkl=parallel "
      "<OBJECTS> "
      "-o <TARGET> "
      "-Xlinker --start-group "  # manually pass to ld.lld
      "<LINK_LIBRARIES> "
      "-Xlinker --end-group"
    )
  elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    string(CONCAT LINK_EXECUTABLE_STR
      "<CMAKE_Fortran_COMPILER> "
      "-nologo "
      "-ipo "  # with -fuse-ld=lld
      "-qopenmp "
      "-qmkl=parallel "
      "<OBJECTS> "
      "-o <TARGET> "
      "-Xlinker --start-group "  # manually pass to ld.lld
      "<LINK_LIBRARIES> "
      "-Xlinker --end-group"
    )
  endif()
endif()

set(CMAKE_Fortran_COMPILE_OBJECT "${COMPILE_OBJECT_STR}")
set(CMAKE_Fortran_LINK_EXECUTABLE "${LINK_EXECUTABLE_STR}")


#-------------------------------------------------------------------------------
# Post-build steps
#-------------------------------------------------------------------------------
# set stack size
if(WIN32)
  add_custom_command(TARGET TRESC POST_BUILD
    COMMAND editbin /STACK:33554432 $<TARGET_FILE:TRESC>
    COMMENT "Setting stack size"
  )
endif()