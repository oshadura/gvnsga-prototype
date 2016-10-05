#pragma once

#ifndef __CMSGEANTV__
#define __CMSGEANTV__

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
//#include "Memory.h"

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

#ifndef COPROCESSOR_REQUEST
#define COPROCESSOR_REQUEST false
#endif

using namespace Geant;

namespace geantvmoop {

class CMSGeantV : public Functions<CMSGeantV> {

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
    static bool monitor = false, score = false, debug = false;
    const char *cms("cms2015.root");
    const char  *xsec("xsec_FTFP_BERT_G496p02_1mev.root");
    const char *fstate("fstate_FTFP_BERT_G496p02_1mev.root");
    std::string  hepmc_event_filename("pp14TeVminbias.root");
    bool performance = true;
    bool coprocessor = COPROCESSOR_REQUEST;
    int n_threads = fParameters[0];
   #ifdef USE_ROOT
  TGeoManager::Import(cms);
#else

#endif

  TaskBroker *broker = nullptr;
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
  
  config->fNtotal = fParameters[0];
  config->fNbuff = fParameters[1];
  // Default value is 1. (0.1 Tesla)
  config->fBmag = 40.; // 4 Tesla

  // Enable use of RK integration in field for charged particles
  config->fUseRungeKutta = false;
  // prop->fEpsilonRK = 0.001;  // Revised / reduced accuracy - vs. 0.0003 default

  config->fNminThreshold=5 * n_threads;
  config->fUseMonitoring = false;
  config->fNaverage = 500;

  config->SetMonitored(GeantConfig::kMonQueue, monitor);
  config->SetMonitored(GeantConfig::kMonMemory, monitor);
  config->SetMonitored(GeantConfig::kMonBasketsPerVol, monitor);
  config->SetMonitored(GeantConfig::kMonVectors, monitor);
  config->SetMonitored(GeantConfig::kMonConcurrency, monitor);
  config->SetMonitored(GeantConfig::kMonTracksPerEvent, monitor);
  // Threshold for prioritizing events (tunable [0, 1], normally <0.1)
  // If set to 0 takes the default value of 0.01
  config->fPriorityThr = fParameters[2]/100;

  // Initial vector size, this is no longer an important model parameter,
  // because is gets dynamically modified to accomodate the track flow
  config->fNperBasket = fParameters[3]*32; // Initial vector size

  // This is now the most important parameter for memory considerations
  config->fMaxPerBasket = fParameters[4]*32;

  // Maximum user memory limit [MB]
  config->fMaxRes = 4000;
  if (config) config->fMaxRes = 0;
  config->fEmin = 0.001; // [1 MeV] energy cut
  config->fEmax = 0.01;  // 10 MeV
  if (debug) {
    config->fUseDebug = true;
    config->fDebugTrk = 1;
    //propagator->fDebugEvt = 0;
    //propagator->fDebugStp = 0;
    //propagator->fDebugRep = 10;
  }
  config->fUseMonitoring = false;

  // Set threshold for tracks to be reused in the same volume
  config->fNminReuse = fParameters[5]*10000;

  // Activate standard scoring   
  config->fUseStdScoring = true;
  if (performance) config->fUseStdScoring = false;
  config->fLearnSteps = fParameters[6]*1000;
  if (performance) config->fLearnSteps = 0;

  // Activate I/O
  config->fFillTree = false;
  config->fTreeSizeWriteThreshold = 100000;
  // Activate old version of single thread serialization/reading
  //   config->fConcurrentWrite = false;

  // Create run manager
  GeantRunManager *runMgr = new GeantRunManager(1, n_threads, config);
  if (broker) runMgr->SetCoprocessorBroker(broker);
  // Create the tab. phys process.
  runMgr->SetPhysicsProcess( new TTabPhysProcess("tab_phys", xsec, fstate));

#ifdef USE_VECGEOM_NAVIGATOR
#ifdef USE_ROOT
//  runMgr->LoadVecGeomGeometry();
#else
//  runMgr->LoadGeometry(cms_geometry_filename.c_str());
#endif
#endif

  if (hepmc_event_filename.empty()) {
    runMgr->SetPrimaryGenerator( new GunGenerator(config->fNaverage, 11, config->fEmax, -8, 0, 0, 1, 0, 0) );
  } else {
    // propagator->fPrimaryGenerator->SetEtaRange(-2.,2.);
    // propagator->fPrimaryGenerator->SetMomRange(0.,0.5);
    // propagator->fPrimaryGenerator = new HepMCGenerator("pp14TeVminbias.hepmc3");
    runMgr->SetPrimaryGenerator( new HepMCGenerator(hepmc_event_filename) );
  }

  CMSApplication *CMSApp = new CMSApplication(runMgr);
  runMgr->SetUserApplication( CMSApp );
  if (score) {
    CMSApp->SetScoreType(CMSApplication::kScore);
  } else {
    CMSApp->SetScoreType(CMSApplication::kNoScore);
  }
#ifdef GEANT_TBB
  if (tbbmode)
    runMgr->SetTaskMgr( new TaskMgrTBB() );
#endif

  runMgr->RunSimulation();
    fFitness.push_back(timer.Elapsed());
    //fFitness.push_back(-(runMgr->fNprimaries.load()/timer.Elapsed()));
#ifdef ENABLE_PERF
    ////size_t peakSize = getPeakRSS();
    //fFitness.push_back(perfcontrol.GetNICS());
    //fFitness.push_back(perfcontrol.GetNCS());
    //fFitness.push_back(perfcontrol.GetNC());
    //fFitness.push_back(perfcontrol.GetNI());
    //fFitness.push_back(perfcontrol.GetNBM());
    //fFitness.push_back(perfcontrol.GetNDC());
    //fFitness.push_back(perfcontrol.GetNIC());
    //fFitness.push_back(perfcontrol.GetNB());
    ////fFitness.push_back(peakSize);
    //fFitness.push_back(cpumgr.GetCurrentValueCPU());
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
    for (int i = 0; i < 7; ++i)
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

