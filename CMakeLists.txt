project(CommonAlgo)
#cmake_minimum_required(VERSION 3.13)

find_package(EventAnalysis)
find_package(HerdSoftware)
find_package(ROOT)
include_directories(${ROOT_INCLUDE_DIRS})


add_library(CommonAlgo  SHARED 
			MCtruthProcess.cpp
                        CaloGeomFidVolume.cpp
           )

target_link_libraries(CommonAlgo EACore EAData EAAlgorithm EAUtils EAAnalysis HerdDataObjects ${ROOT_LIBRARIES})

