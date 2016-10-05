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

#ifndef COPROCESSOR_REQUEST
#define COPROCESSOR_REQUEST false
#endif

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
    const char *cms_geometry_filename("cms2015.root");
    const char  *xsec_filename("xsec_FTFP_BERT_G496p02_1mev.root");
    const char *fstate_filename("fstate_FTFP_BERT_G496p02_1mev.root");
    std::string  hepmc_event_filename("pp14TeVminbias.root");
    bool performance = true;
    bool coprocessor = COPROCESSOR_REQUEST;
    int nthreads = fParameters[0];
    printf("Debugging Run.C: thread value = %d\n", nthreads);
    int ntotal = fParameters[1];
    printf("Debugging Run.C: all events value = %d\n", ntotal);
    int nbuffered = fParameters[2];
    printf("Debugging Run.C: buffered particles value = %d\n", nbuffered);
    TGeoManager::Import(cms_geometry_filename);
    TaskBroker *broker = nullptr;
    if (coprocessor) {
#ifdef GEANTCUDA_REPLACE
      CoprocessorBroker *gpuBroker = new CoprocessorBroker();
      gpuBroker->CudaSetup(32, 128, 1);
      broker = gpuBroker;
      nthreads += gpuBroker->GetNstream() + 1;
#else
      std::cerr
          << "Error: Coprocessor processing requested but support was not "
             "enabled\n";
#endif
    }
    WorkloadManager *wmanager = WorkloadManager::Instance(ntotal);
    std::cout
        << "-=======================GeantPropagator=======================-"
        << std::endl;
    GeantPropagator *prop =
        GeantPropagator::Instance(nbuffered, nthreads);
    if (broker)
      prop->SetTaskBroker(broker);
    // Monitor different features
    wmanager->SetNminThreshold(5 * nthreads);
    wmanager->SetMonitored(GeantPropagator::kMonQueue, true & (!performance));
    wmanager->SetMonitored(GeantPropagator::kMonMemory, false & (!performance));
    wmanager->SetMonitored(GeantPropagator::kMonBasketsPerVol,
                       false & (!performance));
    wmanager->SetMonitored(GeantPropagator::kMonVectors, false & (!performance));
    wmanager->SetMonitored(GeantPropagator::kMonConcurrency,
                       false & (!performance));
    wmanager->SetMonitored(GeantPropagator::kMonTracksPerEvent,
                       false & (!performance));
    bool graphics = (prop->GetMonFeatures()) ? true : false;
    prop->fUseMonitoring = graphics;
    prop->fNaverage = 500; // Average number of tracks per event
    // Threshold for prioritizing events (tunable [0, 1], normally <0.1)
    // If set to 0 takes the default value of 0.01
    prop->fPriorityThr = fParameters[3];
    printf("Debugging RunCMS.C: priority value = %f\n", prop->fPriorityThr);
    // Initial vector size, this is no longer an important model parameter,
    // because is gets dynamically modified to accomodate the track flow
    prop->fNperBasket = fParameters[4];
    printf("Debugging RunCMS.C: vector value = %d\n", prop->fNperBasket);
    // This is now the most important parameter for memory considerations
    prop->fMaxPerBasket = 256; // Maximum vector size (tunable)
    int max_memory = 4000;
    prop->fMaxRes = max_memory;
    if (performance) prop->fMaxRes = 0;
    prop->fEmin = 0.001;       // [3 KeV] energy cut
    prop->fEmax = 0.01; // [30MeV] used for now to select particle gun energy
    // Create the tab. phys process.
    prop->LoadGeometry(cms_geometry_filename);
    std::cout
        << "-=======================TTabPhysProcess=======================-"
        << std::endl;
    prop->fProcess = new TTabPhysProcess("tab_phys", xsec_filename, fstate_filename);
    // for vector physics -OFF now
    // prop->fVectorPhysicsProcess = new GVectorPhysicsProcess(prop->fEmin,
    // nthreads);
    std::cout << "-=======================GunGenerator=======================-"
              << std::endl;
    //prop->fPrimaryGenerator =
    //    new GunGenerator(prop->fNaverage, 11, prop->fEmax, -8, 0, 0, 1, 0, 0);
    prop->fPrimaryGenerator = new HepMCGenerator(hepmc_event_filename);
    // Number of steps for learning phase (tunable [0, 1e6])
    // if set to 0 disable learning phase
    prop->fLearnSteps = fParameters[5];
    printf("Debugging Run.C: learning steps value = %d\n", prop->fLearnSteps);
    if (performance)
      prop->fLearnSteps = 0;
    std::cout
        << "-=======================ExN03Application=======================-"
        << std::endl;
    prop->fFillTree = false;
    prop->fTreeSizeWriteThreshold = 100000;
    //prop->fApplication = new ExN03Application();
    CMSApplication *CMSApp = new CMSApplication();
    if (score) {
      CMSApp->SetScoreType(CMSApplication::kScore);
    } else {
     CMSApp->SetScoreType(CMSApplication::kNoScore);
    }
    prop->fApplication = CMSApp;
    // Activate I/O
    prop->fFillTree = false;
    // Activate old version of single thread serialization/reading
    // prop->fConcurrentWrite = false;
    // Activate debugging using -DBUG_HUNT=ON in your cmake build
    prop->fDebugEvt = 0;
    prop->fDebugTrk = 0;
    prop->fDebugStp = 0;
    prop->fDebugRep = 10;
    // Activate standard scoring
    prop->fUseStdScoring = true;
    if (performance)
      prop->fUseStdScoring = false;
    // Monitor the application
    //prop->fUseAppMonitoring = false;
    prop->PropagatorGeom(cms_geometry_filename, nthreads, graphics);
#ifdef ENABLE_PERF
    perfcontrol.Stop();
    timer.Stop();
#endif
    fFitness.push_back(timer.Elapsed());
    fFitness.push_back(-(prop->fNprimaries.load()/timer.Elapsed()));
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
    //delete prop;
    std::cout << "Vector output for evaluation function: ";
    for (auto i : fFitness)
      std::cout << i << ' ';
    std::cout << ' ' << std::endl;
    return fFitness;
  }

  static Input GetInput() {
    Input vector;
    for (int i = 0; i < 6; ++i)
      vector.push_back(GADouble(3, 12));
    return vector;
  }
#ifndef ENABLE_PERF
  static Output GetOutput() { return std::vector<double>(2); }
#else
  static Output GetOutput() { return std::vector<double>(4); }
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

