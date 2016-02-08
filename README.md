# gvnsga-prototype project

First step using hardcoded version with NSGA2 as genetic algorithm for small number of objectited functions (RT & memory)

Deb, K., Pratap. A, Agarwal, S., and Meyarivan, T. (2002). A fast and elitist multiobjective genetic algorithm: NSGA-II. IEEE Transaction on Evolutionary Computation, 6(2), 181-197. 

2.12.15
Simplified version with hardcoded runCMS.C function and Genes/Fitness vectors

What is needed to be run:

-> GeantV configuration cmake file (changes in GeantV cmake file)

How to run:

-> geantbuild.cern.ch

cmake ../ -DGeantV_DIR=/home/geant/Install/GeantV/cmake/ -DVc_DIR=/home/geant/Install/VC1.0/lib/cmake/Vc/ -DHEPMC_ROOT_DIR=/afs/cern.ch/work/g/geant/jenkins/workspace/olwork21-gcc491-geant/HepMC/ -DCMAKE_CXX_FLAGS="-fabi-version=6 -fdiagnostics-color"

####### New version #######

cmake ../ -DCMAKE_INSTALL_PREFIX=/home/geant/Install/gvnsga-prototype -DVc_DIR=/home/geant/Install/Vc/lib/cmake/Vc/ -DGeantV_DIR=/home/geant/Install/GeantV/cmake/ -DCMAKE_PREFIX_PATH=/home/geant/Install/HepMC/ -DCMAKE_CXX_FLAGS="-fdiagnostics-color=auto -fpermissive -fabi-version=6 -fext-numeric-literals" -DVecGeom_DIR=/home/geant/Install/VecGeom/lib/CMake/ -DCMAKE_BUILD_TYPE=Debug -DBOOST_ROOT=/home/geant/Install/boost/

############ BUILD with GeantV integration ##############

cmake ../ -DCMAKE_INSTALL_PREFIX=/home/geant/Install/gvnsga-prototype -DVc_DIR=/home/geant/Install/Vc/lib/cmake/Vc/ -DGeantV_DIR=/home/geant/Install/GeantV/cmake/ -DCMAKE_PREFIX_PATH=/home/geant/Install/HepMC/ -DCMAKE_CXX_FLAGS="-fdiagnostics-color=auto -fpermissive -fabi-version=6 -fext-numeric-literals" -DVecGeom_DIR=/home/geant/Install/VecGeom/lib/CMake/VecGeom/ -DCMAKE_BUILD_TYPE=Debug -DBOOST_ROOT=/home/geant/Install/boost/ -DENABLE_GEANTV=ON -DNUMERIC_LIB=OFF

############## BUILD with numerical examples #############

cmake ../ -DCMAKE_INSTALL_PREFIX=/home/geant/Install/gvnsga-prototype -DVc_DIR=/home/geant/Install/Vc/lib/cmake/Vc/ -DGeantV_DIR=/home/geant/Install/GeantV/cmake/ -DCMAKE_PREFIX_PATH=/home/geant/Install/HepMC/ -DCMAKE_CXX_FLAGS="-fdiagnostics-color=auto -fpermissive -fabi-version=6 -fext-numeric-literals" -DVecGeom_DIR=/home/geant/Install/VecGeom/lib/CMake/VecGeom/ -DCMAKE_BUILD_TYPE=Debug -DBOOST_ROOT=/home/geant/Install/boost/ -DENABLE_GEANTV=OFF -DNUMERIC_LIB=ON
