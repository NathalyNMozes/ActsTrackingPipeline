cmake_minimum_required(VERSION 3.14)

# LCG sets CPATH, which gets treated like -I by the compiler. We want to ignore
# warnings from libraries, by unsetting it here, it gets processed by the usual
# target_include_directories call, resulting in the desired -isystem flag.
unset(ENV{CPATH})

project(TrackingPipeline LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

# set the install directories
# to be placed in the build directory
set(CMAKE_INSTALL_BINDIR bin)
set(CMAKE_INSTALL_LIBDIR lib)

# place build products in `<build>/bin` and `<build>/lib` for simple use
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}")

# TBB, ROOT::Core are required for the
# algorithm execution
find_package(
    TBB 
    COMPONENTS 
    tbb 
    REQUIRED
)
find_package(
    ROOT 
    COMPONENTS 
    Core
    RIO
    Tree
    Hist
    HistPainter
    Gpad
    Net
    Physics
    GeomPainter
    REQUIRED
)

find_package(
  OpenMP 
  COMPONENTS
  CXX
  REQUIRED
)

# EUDAQ configuration (package or local build)
# First: try to use a proper EUDAQ CMake package 
# find_package(
#   eudaq
#   2.5.2
#   COMPONENTS
#     stave
#     detector_event
#   QUIET
# )

# set(HAVE_EUDAQ_PACKAGE OFF)
# if (eudaq_FOUND)
#   set(HAVE_EUDAQ_PACKAGE ON)
#   message(STATUS "Found EUDAQ package at ${eudaq_DIR}")
#   include("${eudaq_DIR}/eudaqStaveConfig.cmake")
# endif()

# # If no full package, allow a local minimal build with dict+headers (from_eudaq)
# if (NOT HAVE_EUDAQ_PACKAGE)
#   message(STATUS
#     "EUDAQ package not found, trying local minimal EUDAQ build via EUDAQ_ROOT")

#   set(EUDAQ_ROOT "" CACHE PATH
#       "Path to a minimal EUDAQ build (containing inc/ and build/libdet_event_dict.*)")

#   if (NOT EUDAQ_ROOT)
#     message(FATAL_ERROR
#       "EUDAQ not found.\n"
#       "Either install EUDAQ and set eudaq_DIR, or\n"
#       "set EUDAQ_ROOT to your local minimal build (e.g. /path/to/from_eudaq).")
#   endif()

#   set(EUDAQ_INC_DIR "${EUDAQ_ROOT}/inc" CACHE PATH "EUDAQ include directory")
#   set(EUDAQ_DICT    "${EUDAQ_ROOT}/build/libdet_event_dict.dylib"
#                     CACHE FILEPATH "EUDAQ detector_event dictionary library")

#   if (NOT EXISTS "${EUDAQ_INC_DIR}")
#     message(FATAL_ERROR "EUDAQ_INC_DIR does not exist: ${EUDAQ_INC_DIR}")
#   endif()
#   if (NOT EXISTS "${EUDAQ_DICT}")
#     message(FATAL_ERROR "EUDAQ_DICT does not exist: ${EUDAQ_DICT}")
#   endif()
# endif()

# Dictionary for EUDAQ detector_event
# add_library(DetEventDict SHARED)

# target_include_directories(
#     DetEventDict
#     PUBLIC
#     $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
# )

# ROOT_GENERATE_DICTIONARY(
#     G__DetEventDict
#     ./include/TrackingPipeline/from_eudaq/DetectorEvent.h
#     MODULE
#     DetEventDict
#     LINKDEF
#     ./include/TrackingPipeline/from_eudaq/LinkDef.h
# )

# target_link_libraries(
#     DetEventDict
#     PRIVATE
#     ROOT::RIO
#     ROOT::Physics
#     ROOT::GeomPainter
# )


# Dictionary for apollonIo detector_event
add_library(DetEventClusteredDict SHARED)
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include/TrackingPipeline/Preprocessing
)
ROOT_GENERATE_DICTIONARY(
    G__DetEventClusteredDict 
    ./include/TrackingPipeline/Preprocessing/DetectorEvent.hpp
    MODULE 
    DetEventClusteredDict 
    LINKDEF 
    ./include/TrackingPipeline/Preprocessing/LinkDef.hpp
)
target_link_libraries(
    DetEventClusteredDict
    PRIVATE 
    ROOT::RIO
    ROOT::Physics
    ROOT::GeomPainter
)

# Dictionary for the sim data
add_library(SimEventDict SHARED)
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include/TrackingPipeline/Io/detail
)
ROOT_GENERATE_DICTIONARY(
    G__SimEventDict
    ./include/TrackingPipeline/Io/detail/ROOTDictDefs.h
    MODULE 
    SimEventDict 
    LINKDEF 
    ./include/TrackingPipeline/Io/detail/LinkDef.h
)
target_link_libraries(
    SimEventDict 
    PRIVATE 
    ROOT::RIO
    ROOT::Physics
    ROOT::GeomPainter
    DetEventClusteredDict
    # eudaq_preprocessor::detector_event
)

# Geant4 is required for the geometry
find_package(Geant4 REQUIRED)

# Boost is required for the Acts
find_package(Boost REQUIRED)

# The core components and the Algo
# infrastructure are required
find_package(
    Acts
    COMPONENTS
    Core
    Fatras
    PluginGeant4
    PluginJson
    PluginFpeMonitoring
    Alignment
    REQUIRED
)

# toml11: header-only TOML library (installed via system package manager)
find_package(toml11 REQUIRED)
# source files will be added later
add_library(
    TrackingPipelineCore SHARED "")

target_link_libraries(
    TrackingPipelineCore
    PUBLIC
    TBB::tbb
    OpenMP::OpenMP_CXX
    ROOT::Core
    ROOT::RIO
    ROOT::Tree
    ROOT::Physics
    ROOT::Hist
    ROOT::Gpad
    ROOT::Net
    G__SimEventDict
    ActsCore
    ActsFatras
    ActsPluginGeant4
    ActsPluginJson
    ActsPluginFpeMonitoring
    ActsAlignment
    DetEventClusteredDict
    # DetEventDict
    toml11::toml11
    ActsPreprocessing
)

target_include_directories(
    TrackingPipelineCore
    PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

# collect the library dependencies
add_subdirectory(src/Geometry)
add_subdirectory(src/Infrastructure)
add_subdirectory(src/Io)
add_subdirectory(src/MagneticField)
add_subdirectory(src/Material)
add_subdirectory(src/Simulation)
add_subdirectory(src/TrackFinding)
add_subdirectory(src/TrackFitting)
add_subdirectory(src/Alignment)
add_subdirectory(src/Run)
add_subdirectory(src/Analysis)
add_subdirectory(src/Preprocessing)


# macro to add a run executables
macro(add_run _name _source)
    # automatically prefix the target name
    set(_target "Run${_name}")
    add_executable(${_target} ${_source})
    # core link is public, so no need 
    # to add the includes
    target_link_libraries(
        ${_target}
        PRIVATE
        TrackingPipelineCore
    )
endmacro()

# run has to be last, so that it can be 
# linked against the library
add_subdirectory(Run)
