project(acceptanceAlgo)
cmake_minimum_required(VERSION 3.13)

find_package(EventAnalysis)
find_package(HerdSoftware)
find_package(ROOT)
include_directories(${ROOT_INCLUDE_DIRS})


add_library(acceptanceAlgo  SHARED 
			GeomAcceptance/MCtruthProcess.cpp
                        GeomAcceptance/CaloGeomFidVolume.cpp
			Calo/CaloGlob.cpp
                        Calo/CaloAxis.cpp
                        Calo/CaloAxisInfo.cpp
                        Calo/CaloTest.cpp
                        Histo/mcEnergyHisto.cpp
                        Histo/mcGenSpectrum.cpp
                        Histo/mcEnergyHisto3D.cpp
           )

target_link_libraries(acceptanceAlgo EACore EAData EAAlgorithm EAUtils EAAnalysis HerdDataObjects ${ROOT_LIBRARIES})
