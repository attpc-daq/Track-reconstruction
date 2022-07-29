
#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(HoughTransform)

#----------------------------------------------------------------------------
# Cmake Path for installation
#
#CMAKE_INSTALL_PREFIX
#----------------------------------------------------------------------------
# Load some basic macros which are needed later on
include(../FindROOT.cmake)

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
# Setup ROOT include directories

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${ROOT_INCLUDE_DIR})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hpp)

#----------------------------------------------------------------------------
# Set everything needed for the root dictonary and create the
# dictionary
set(TrackRec_HEADERS  ${PROJECT_SOURCE_DIR}/include/HoughTransform2D.hpp)

set(TrackRec_LINKDEF ${PROJECT_SOURCE_DIR}/src/FitterLinkDef.hh )
set(TrackRec_DICTIONARY FitterDict.cc)
ROOT_GENERATE_DICTIONARY("${TrackRec_HEADERS}" "${TrackRec_LINKDEF}" "${TrackRec_DICTIONARY}" "${ROOT_INCLUDE_DIR}")

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 and ROOT libraries
#
add_executable(HoughTransform ${PROJECT_SOURCE_DIR}/main.cpp ${TrackRec_DICTIONARY} ${sources} ${headers})
target_link_libraries(HoughTransform  ${ROOT_LIBRARIES}
Spectrum
Geom
EG
Core
Gpad
Graf
Graf3d
Gui
Hist
MathCore
Matrix
Minuit
Net
Physics
Postscript
Rint
RIO
Thread
Tree
TreePlayer)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
# install(TARGETS  Analysis_Angle DESTINATION bin)
install(TARGETS  HoughTransform DESTINATION bin)

