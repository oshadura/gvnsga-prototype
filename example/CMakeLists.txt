#----------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)
#----------------------------------------------------------------------------
project(libGAExamples)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/../cmake/modules")
#----------------------------------------------------------------------------
find_package(ROOT REQUIRED)
find_package(GeantV REQUIRED)
find_package(VecGeom REQUIRED)
#find_package(Vc REQUIRED)
find_package(HepMC REQUIRED)
find_package(Boost)
#---------------------------------------------------------------------------
include_directories(${Vc_INCLUDE_DIR} ${HEPMC_INCLUDE_DIR} ${VECGEOM_INCLUDE_DIR} ${GEANTV_INCLUDE_DIR})
include_directories(${CMAKE_SOURCE_DIR}/libga/inc)
include(${ROOT_USE_FILE})

#---------------------------------------------------------------------------

if(Boost_FOUND)
	include_directories(${Boost_INCLUDE_DIRS})
endif()

#----------------------------------------------------------------------------

message(STATUS "GEANTV_INCLUDE_DIR is ${GEANTV_INCLUDE_DIR}")
message(STATUS "GEANTV_LIBRARIES is ${GEANTV_LIBRARIES}")
message(STATUS "VECGEOM_LIBRARIES is ${VECGEOM_LIBRARIES}")
message(STATUS "VECGEOM_INCLUDE_DIR is ${VECGEOM_INCLUDE_DIR}")

#----------------------------------------------------------------------------
macro(build_executables EXECUTABLES)
  foreach(EXECUTABLE ${EXECUTABLES})
    get_filename_component(TARGET_NAME ${EXECUTABLE} NAME_WE)
    add_executable(${TARGET_NAME} ${EXECUTABLE})
    if(Boost_FOUND)
    	target_link_libraries(${TARGET_NAME} GA ${ROOT_LIBRARIES} Hist RIO Tree Gui ${GEANTV_LIBRARIES} ${VECGEOM_LIBRARIES} ${VC_LIBRARIES} ${HEPMC_LIBRARIES} ${Boost_LIBRARY_DIRS} -lpfm)
    else()
    	target_link_libraries(${TARGET_NAME} GA ${ROOT_LIBRARIES} Hist RIO Tree Gui ${GEANTV_LIBRARIES} ${VECGEOM_LIBRARIES} ${VC_LIBRARIES} ${HEPMC_LIBRARIES} -lpfm)
    endif()
    add_dependencies(${TARGET_NAME} GA)
  endforeach()
endmacro()

#----------------------------------------------------------------------------
set(LIBGA_EXAMPLES
   ${CMAKE_SOURCE_DIR}/example/PopulationTreeDTLZ1.cxx
   ${CMAKE_SOURCE_DIR}/example/PopulationTreeDTLZ2.cxx
   ${CMAKE_SOURCE_DIR}/example/PopulationTreeDTLZ3.cxx
   ${CMAKE_SOURCE_DIR}/example/PopulationTreeRun.cxx
   ${CMAKE_SOURCE_DIR}/example/PopulationTreeRunCMS.cxx
)

#----------------------------------------------------------------------------

#add_executable(PopulationTree ${CMAKE_SOURCE_DIR}/example/PopulationTree.cxx)
#target_link_libraries(PopulationTree GA ${ROOT_LIBRARIES} Hist RIO Tree Gui ${VECGEOM_LIBRARIES} ${HEPMC3_LIBRARIES} ${GEANTV_LIBRARIES})

build_executables("${LIBGA_EXAMPLES}")
install(TARGETS PopulationTreeDTLZ1 PopulationTreeDTLZ2 PopulationTreeDTLZ3 PopulationTreeRun PopulationTreeRunCMS DESTINATION ${CMAKE_BINARY_DIR})
