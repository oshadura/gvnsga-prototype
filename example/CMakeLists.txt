#----------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)
#----------------------------------------------------------------------------
project(libGAExamples)
#list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/../cmake/modules")
#----------------------------------------------------------------------------
#find_package(ROOT REQUIRED)
#find_package(GeantV REQUIRED)
#find_package(VecGeom REQUIRED)
#find_package(Vc REQUIRED)

#----------------------------------------------------------------------------
#if(ENABLE_PERFMON)
#  find_package(PerfMon)
#endif()
#----------------------------------------------------------------------------

if(ENABLE_GEANTV)
  find_package(HepMC)
endif()
#---------------------------------------------------------------------------
include_directories(${Vc_INCLUDE_DIR} ${HEPMC_INCLUDE_DIR} ${VECGEOM_INCLUDE_DIR} ${GEANTV_INCLUDE_DIR} ${PerfMon_INCLUDE_DIR})
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
message(STATUS "PerfMon_LIBRARY is ${PerfMon_LIBRARY}")
message(STATUS "PerfMon_INCLUDE_DIR is ${PerfMon_INCLUDE_DIR}")

#----------------------------------------------------------------------------
macro(build_executables EXECUTABLES)
  foreach(EXECUTABLE ${EXECUTABLES})
    get_filename_component(TARGET_NAME ${EXECUTABLE} NAME_WE)
    add_executable(${TARGET_NAME} ${EXECUTABLE})
    #if(PERFMON_FOUND)
    	target_link_libraries(${TARGET_NAME} GA ${ROOT_LIBRARIES} Hist RIO Tree Gui ${GEANTV_LIBRARIES} ${VECGEOM_LIBRARIES} ${VC_LIBRARIES} ${HEPMC_LIBRARIES} ${Boost_LIBRARY_DIRS} ${PerfMon_LIBRARY} ${MLPACK_LIBRARY})
    #else()
    #	target_link_libraries(${TARGET_NAME} GA ${ROOT_LIBRARIES} Hist RIO Tree Gui ${GEANTV_LIBRARIES} ${VECGEOM_LIBRARIES} ${VC_LIBRARIES} ${HEPMC_LIBRARIES} ${Boost_LIBRARY_DIRS})
    #endif()
    add_dependencies(${TARGET_NAME} GA)
  endforeach()
endmacro()

#----------------------------------------------------------------------------
set(LIBGA_EXAMPLES
   ${CMAKE_SOURCE_DIR}/example/PopulationTreeDTLZ1.cxx
   ${CMAKE_SOURCE_DIR}/example/PopulationTreeDTLZ2.cxx
   ${CMAKE_SOURCE_DIR}/example/PopulationTreeDTLZ3.cxx
   ${CMAKE_SOURCE_DIR}/example/PopulationTreeDTLZ4.cxx
   ${CMAKE_SOURCE_DIR}/example/PopulationTreeDTLZ5.cxx
   ${CMAKE_SOURCE_DIR}/example/PopulationTreeDTLZ6.cxx
   ${CMAKE_SOURCE_DIR}/example/PopulationTreeDTLZ7.cxx
   ${CMAKE_SOURCE_DIR}/example/PopulationTreeRunNSGA2.cxx
   ${CMAKE_SOURCE_DIR}/example/PopulationTreeRunCMS.cxx
   ${CMAKE_SOURCE_DIR}/example/PopulationTreeRunCMAES.cxx
   #${CMAKE_SOURCE_DIR}/example/PopulationCMSTest.cxx
)

#----------------------------------------------------------------------------

#add_executable(PopulationTree ${CMAKE_SOURCE_DIR}/example/PopulationTree.cxx)
#target_link_libraries(PopulationTree GA ${ROOT_LIBRARIES} Hist RIO Tree Gui ${VECGEOM_LIBRARIES} ${HEPMC3_LIBRARIES} ${GEANTV_LIBRARIES})

build_executables("${LIBGA_EXAMPLES}")
install(TARGETS
  PopulationTreeDTLZ1
  PopulationTreeDTLZ2
  PopulationTreeDTLZ3
  PopulationTreeDTLZ4
  PopulationTreeDTLZ5
  PopulationTreeDTLZ6
  PopulationTreeDTLZ7
  PopulationTreeRunNSGA2
  PopulationTreeRunCMS
  #PopulationCMSTest
 DESTINATION ${CMAKE_BINARY_DIR})
