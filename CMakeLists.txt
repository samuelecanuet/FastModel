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
