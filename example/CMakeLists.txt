#----------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)
#----------------------------------------------------------------------------
project(PopulationTree)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/../cmake/modules")

find_package(ROOT REQUIRED)
find_package(GeantV REQUIRED)
find_package(VecGeom REQUIRED)
find_package(Vc REQUIRED)
find_package(HepMC REQUIRED)
include_directories(${Vc_INCLUDE_DIR} ${HEPMC_INCLUDE_DIR} ${VECGEOM_INCLUDE_DIR} ${GEANTV_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/libga/inc)
include(${ROOT_USE_FILE})
#----------------------------------------------------------------------------

message(STATUS "GEANTV_INCLUDE_DIR is -${GEANTV_INCLUDE_DIR}-")
message(STATUS "GEANTV_LIBRARIES is -${GEANTV_LIBRARIES}-")
message(STATUS "VECGEOM_LIBRARIES is -${VECGEOM_LIBRARIES}-")
message(STATUS "VECGEOM_INCLUDE_DIR is -${VECGEOM_INCLUDE_DIR}-")

#----------------------------------------------------------------------------
add_executable(PopulationTree ${CMAKE_SOURCE_DIR}/example/PopulationTree.cxx)
target_link_libraries(PopulationTree GA ${ROOT_LIBRARIES} Hist RIO Tree Gui ${VECGEOM_LIBRARIES} ${HEPMC3_LIBRARIES} ${GEANTV_LIBRARIES})
install(TARGETS PopulationTree DESTINATION ${CMAKE_BINARY_DIR})
