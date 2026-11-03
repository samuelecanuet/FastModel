cmake_minimum_required(VERSION 3.16)

# This directory can be built directly for a compile check or added to a
# parent simulation after it has found Geant4.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  project(UniformFieldFastModel LANGUAGES CXX)
  find_package(Geant4 REQUIRED)
  include(${Geant4_USE_FILE})
endif()

add_library(UniformFieldFastModel
  src/UniformFieldFastGeometry.cc
  src/UniformFieldFastMessenger.cc
  src/UniformFieldFastSimulationModel.cc
  src/UniformFieldFastPropagator.cc
  include/UniformFieldFastGeometry.hh
  include/UniformFieldFastMessenger.hh
  include/UniformFieldFastSimulationModel.hh
  include/UniformFieldFastPropagator.hh
)
add_library(FastModel::UniformField ALIAS UniformFieldFastModel)
target_include_directories(UniformFieldFastModel
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
)
target_link_libraries(UniformFieldFastModel PUBLIC ${Geant4_LIBRARIES})
target_compile_features(UniformFieldFastModel PUBLIC cxx_std_17)
