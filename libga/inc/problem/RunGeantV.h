#pragma once

#ifndef __RUNGEANTV__
#define __RUNGEANTV__

#include <algorithm> // std::copy
#include <cmath>
#include <iostream> // std::cout
#include <iterator> // std::ostream_iterator
#include <utility>
#include <vector> // std::vector

#ifdef ENABLE_GEANTVVV

#include "algorithms/GANSGA2.h"
#include "generic/Functions.h"
#include "generic/Functions.h"
#include "generic/GADouble.h"
#include "generic/GAVector.h"
#include "generic/Population.h"
#include "generic/Population.h"
#include "generic/TGenes.h"
#include "generic/TGenes.h"
#include "instrumentation/GeantVFitness.h"
#include "instrumentation/GeantVFitness.h"
#include "instrumentation/CPUManager.h"
#include "output/HistogramManager.h"
#include "output/HistogramManager.h"
#include <boost/math/constants/constants.hpp>

#include <cmath>
#include <utility>

#ifdef ENABLE_PERF
#include "PFMWatch.h"
#endif

#include "Rtypes.h"
#include "TGeoManager.h"
#include "Memory.h"

#include "CMSApplication.h"
#include "ExN03Application.h"
#include "GeantPropagator.h"
#include "GeantVApplication.h"
#include "GunGenerator.h"
#include "HepMCGenerator.h"
#include "TTabPhysProcess.h"
#include "TaskBroker.h"
#include "WorkloadManager.h"
#include "base/Stopwatch.h"
#include "GeantRunManager.h"

#ifdef GEANT_TBB
#include "TaskMgrTBB.h"
#endif

#ifndef COPROCESSOR_REQUEST
#define COPROCESSOR_REQUEST false
#endif

using namespace Geant;

namespace geantvmoop {

class RunGeantV : public Functions<RunGeantV> {

public:
  typedef GAVector<GADouble> Input;

  typedef std::vector<double> Output;

  /*
  private:
    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive &ar, const unsigned int version) {
      ar & boost::serialization::base_object<Functions<RunGeantV>>(*this);
    }

  public:
  */

  static Output Evaluate(const Input &individual) {
    // Converting values
    std::vector<double> fFitness;
    std::vector<double> fParameters;
    fFitness.reserve(individual.size());
    fParameters.reserve(individual.size());
    for (auto parameter : individual)
      fParameters.push_back(parameter.GetGAValue());
#ifdef ENABLE_PERF
    PFMWatch perfcontrol;
    perfcontrol.Start();
#endif
    CPUManager cpumgr;
    cpumgr.InitCPU();
    vecgeom::Stopwatch timer;
    timer.Start();
    const char *geomfile = "ExN03.root";
    const char *xsec = "xsec_FTFP_BERT.root";
    const char *fstate = "fstate_FTFP_BERT.root";
    bool performance = true;
    bool monitor = false;
    bool debug = false;

    std::cout << "Vector input for evaluation function: ";
    for (auto i : fParameters)
      std::cout << i << ' ';
    std::cout << ' ' << std::endl;
    int nthreads = fParameters[0];
    bool coprocessor = COPROCESSOR_REQUEST;
  Geant::TaskBroker *broker = nullptr;
  if (coprocessor) {
#ifdef GEANTCUDA_REPLACE
    CoprocessorBroker *gpuBroker = new CoprocessorBroker();
    gpuBroker->CudaSetup(32,128,1);
    broker = gpuBroker;
    nthreads += gpuBroker->GetNstream()+1;
#else
    std::cerr << "Error: Coprocessor processing requested but support was not enabled\n";
#endif
  }

  GeantConfig* config=new GeantConfig();

  
//  TGeoManager::Import(exn03_geometry_filename.c_str());
  config->fGeomFileName = geomfile;
  config->fNtotal =fParameters[1];
  config->fNbuff = fParameters[1] - 1;
  config->fUseMonitoring = monitor;
  config->fNminThreshold=5*nthreads;
  config->SetMonitored(GeantConfig::kMonQueue, monitor);
  config->SetMonitored(GeantConfig::kMonMemory, monitor);
  config->SetMonitored(GeantConfig::kMonBasketsPerVol, monitor);
  config->SetMonitored(GeantConfig::kMonVectors, monitor);
  config->SetMonitored(GeantConfig::kMonConcurrency, monitor);
  config->SetMonitored(GeantConfig::kMonTracksPerEvent, monitor);
  config->fNaverage = 500;   // Average number of tracks per event
  
  // Threshold for prioritizing events (tunable [0, 1], normally <0.1)
  // If set to 0 takes the default value of 0.01
  config->fPriorityThr = fParameters[2]/100;

  // Initial vector size, this is no longer an important model parameter,
  // because is gets dynamically modified to accomodate the track flow
  config->fNperBasket = fParameters[3]*32;   // Initial vector size (tunable)

  // This is now the most important parameter for memory considerations
  config->fMaxPerBasket = fParameters[4]*32;   // Maximum vector size (tunable)
  config->fEmin = 3.E-6; // [3 KeV] energy cut
  config->fEmax = 0.03;  // [30MeV] used for now to select particle gun energy

   // Number of steps for learning phase (tunable [0, 1e6])
   // if set to 0 disable learning phase
  config->fLearnSteps = fParameters[5]*1000;
  if (performance) config->fLearnSteps = 0;
   // Activate I/O
  config->fFillTree = false;
   // Activate debugging using -DBUG_HUNT=ON in your cmake build
  if (debug) {
    config->fUseDebug = true;
    config->fDebugTrk = 1;
  }
// Activate standard scoring   
  config->fUseStdScoring = true;
  if (performance) config->fUseStdScoring = false;
  // Monitor the application
  config->fUseAppMonitoring = false;

  // Set threshold for tracks to be reused in the same volume
  config->fNminReuse = fParameters[6]*10000;
  int n_propagators;
  // Create run manager
  if(nthreads == fParameters[7])
   n_propagators = fParameters[7] + 1;
  else
   n_propagators = fParameters[7];
  GeantRunManager *runMgr = new GeantRunManager(n_propagators, nthreads, config);
  if (broker) runMgr->SetCoprocessorBroker(broker);
  // Create the tab. phys process.
  runMgr->SetPhysicsProcess( new TTabPhysProcess("tab_phys", xsec, fstate));
  
// Create the tab. phys process.
#ifdef USE_VECGEOM_NAVIGATOR
//  runMgr->LoadVecGeomGeometry();
#endif

  // for vector physics -OFF now
  // runMgr->SetVectorPhysicsProcess(new GVectorPhysicsProcess(config->fEmin, nthreads));
  runMgr->SetPrimaryGenerator( new GunGenerator(config->fNaverage, 11, config->fEmax, -8, 0, 0, 1, 0, 0) );
  runMgr->SetUserApplication ( new ExN03Application(runMgr) );
#ifdef GEANT_TBB
  if (tbbmode)
    runMgr->SetTaskMgr( new TaskMgrTBB() );
#endif
  
  runMgr->RunSimulation();
//  propagator->PropagatorGeom(exn03_geometry_filename.c_str(), n_threads, monitor);
#ifdef ENABLE_PERF
    perfcontrol.Stop();
    timer.Stop();
#endif
    fFitness.push_back(timer.Elapsed());
    fFitness.push_back(-(runMgr->GetNprimaries()/timer.Elapsed()));
#ifdef ENABLE_PERF
    size_t peakSize = getPeakRSS();
    //fFitness.push_back(perfcontrol.GetNICS());
    //fFitness.push_back(perfcontrol.GetNCS());
    //fFitness.push_back(perfcontrol.GetNC());
    //fFitness.push_back(perfcontrol.GetNI());
    //fFitness.push_back(perfcontrol.GetNBM());
    //fFitness.push_back(perfcontrol.GetNDC());
    //fFitness.push_back(perfcontrol.GetNIC());
    //fFitness.push_back(perfcontrol.GetNB());
    fFitness.push_back(peakSize);
    fFitness.push_back(cpumgr.GetCurrentValueCPU());
    perfcontrol.printSummary();
#endif
    //delete runMgr;
    std::cout << "Vector output for evaluation function: ";
    for (auto i : fFitness)
      std::cout << i << ' ';
    std::cout << ' ' << std::endl;
    return fFitness;
  }

  static Input GetInput() {
    Input vector;
    for (int i = 0; i < 8; ++i)
      vector.push_back(GADouble(3, 12));
    return vector;
  }
#ifndef ENABLE_PERF
  static Output GetOutput() { return std::vector<double>(2); }
#else
  static Output GetOutput() { return std::vector<double>(3); }
#endif

  // ROOT Fitting to true Pareto front
  static Double_t TruePF(Double_t *x, Double_t *parameter) {
    Double_t value =
        parameter[0] * x[0] + parameter[1] * x[1] + parameter[2] * x[2] - 0.5;
    return value;
  }
};
}

#endif
#endif
