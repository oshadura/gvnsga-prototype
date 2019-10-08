# GV-NSGAII prototype project

Based on paper of Deb, K., Pratap. A, Agarwal, S., and Meyarivan, T. (2002). A fast and elitist multiobjective genetic algorithm: NSGA-II. IEEE Transaction on Evolutionary Computation, 6(2), 181-197. 

What is needed to be build before:
* GeantV configuration cmake file (changes in GeantV cmake file)

# How to run for example on geantbuild.cern.ch

`cmake ../ -DGeantV_DIR=/home/geant/Install/GeantV/cmake/ -DVc_DIR=/home/geant/Install/VC1.0/lib/cmake/Vc/ -DHEPMC_ROOT_DIR=/afs/cern.ch/work/g/geant/jenkins/workspace/olwork21-gcc491-geant/HepMC/ -DCMAKE_CXX_FLAGS="-fabi-version=6 -fdiagnostics-color"`

# Standalone version

`cmake ../ -DCMAKE_INSTALL_PREFIX=/home/geant/Install/gvnsga-prototype -DVc_DIR=/home/geant/Install/Vc/lib/cmake/Vc/ -DGeantV_DIR=/home/geant/Install/GeantV/cmake/ -DCMAKE_PREFIX_PATH=/home/geant/Install/HepMC/ -DCMAKE_CXX_FLAGS="-fdiagnostics-color=auto -fpermissive -fabi-version=6 -fext-numeric-literals" -DVecGeom_DIR=/home/geant/Install/VecGeom/lib/CMake/ -DCMAKE_BUILD_TYPE=Debug -DBOOST_ROOT=/home/geant/Install/boost/`

# Build with GeantV integration

`cmake ../ -DCMAKE_INSTALL_PREFIX=/home/geant/Install/gvnsga-prototype -DVc_DIR=/home/geant/Install/Vc/lib/cmake/Vc/ -DGeantV_DIR=/home/geant/Install/GeantV/cmake/ -DCMAKE_PREFIX_PATH=/home/geant/Install/HepMC/ -DCMAKE_CXX_FLAGS="-fdiagnostics-color=auto -fpermissive -fabi-version=6 -fext-numeric-literals" -DVecGeom_DIR=/home/geant/Install/VecGeom/lib/CMake/VecGeom/ -DCMAKE_BUILD_TYPE=Debug -DBOOST_ROOT=/home/geant/Install/boost/ -DENABLE_GEANTV=ON -DNUMERIC_LIB=OFF`

# Build with numerical examples

`cmake ../ -DCMAKE_INSTALL_PREFIX=/home/geant/Install/gvnsga-prototype -DVc_DIR=/home/geant/Install/Vc/lib/cmake/Vc/ -DGeantV_DIR=/home/geant/Install/GeantV/cmake/ -DCMAKE_PREFIX_PATH=/home/geant/Install/HepMC/ -DCMAKE_CXX_FLAGS="-fdiagnostics-color=auto -fpermissive -fabi-version=6 -fext-numeric-literals" -DVecGeom_DIR=/home/geant/Install/VecGeom/lib/CMake/VecGeom/ -DCMAKE_BUILD_TYPE=Debug -DBOOST_ROOT=/home/geant/Install/boost/ -DENABLE_GEANTV=OFF -DNUMERIC_LIB=ON`

# Dependencies:

* GeantV

* Vc

* VecGeom

* Boost

* Google test

* Additional dependencies: dlib, mlpack, tensorflow, R

