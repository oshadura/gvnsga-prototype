#include <cmath>
#include <iostream>  // std::cout
#include <iterator>  // std::ostream_iterator
#include <vector>    // std::vector
#include <algorithm> // std::copy

#ifdef ENABLE_GEANTV

#include "generic/Population.h"
#include "generic/Functions.h"
#include "output/HistogramManager.h"
#include "generic/TGenes.h"
#include "algorithms/AlgorithmNSGA.h"
#include "instrumentation/GeantVFitness.h"

#ifdef ENABLE_PERFMON
#include "PFMWatch.h"
#endif

#include "Rtypes.h"
#include "TGeoManager.h"

#include "GunGenerator.h"
#include "HepMCGenerator.h"
#include "TaskBroker.h"
#include "WorkloadManager.h"
#include "GeantPropagator.h"
#include "TTabPhysProcess.h"
#include "CMSApplication.h"
#include "ExN03Application.h"
#include "GeantVApplication.h"

#ifndef COPROCESSOR_REQUEST
#define COPROCESSOR_REQUEST false
#endif

// Forget about constrains now!
std::vector<Double_t> CMSApp(Genes<Double_t> &individual) {
  gROOT->Reset();
  GeantVFitness *fitness = new GeantVFitness();
#ifdef ENABLE_PERFMON
  PFMWatch perfcontrol;
  perfcontrol.Start();
#endif
  std::cout << "Lets pass it to GeantV propagator.." << std::endl;
  bool performance = true;
  const char *geomfile = "cms2015.root";
  const char *xsec = "xsec_FTFP_BERT_G496p02_1mev.root";
  const char *fstate = "fstate_FTFP_BERT_G496p02_1mev.root";
  bool coprocessor = COPROCESSOR_REQUEST;
  // int nthreads = ncputhreads;
  int nthreads = individual.GetThread(individual);
  printf("Debugging RunCMS.C: thread value = %d\n", nthreads);
  // Value from individual vector
  int ntotal =
      individual.GetAllev(individual); // Number of events to be transported
  printf("Debugging RunCMS.C: all events value = %d\n", ntotal);
  // Value from individual vector
  int nbuffered = individual.GetBuffev(
      individual); // Number of buffered events (tunable [1,ntotal])
  printf("Debugging RunCMS.C: buffered particles value = %d\n", nbuffered);
  TGeoManager::Import(geomfile);
  TaskBroker *broker = nullptr;
  if (coprocessor) {
#ifdef GEANTCUDA_REPLACE
    CoprocessorBroker *gpuBroker = new CoprocessorBroker();
    gpuBroker->CudaSetup(32, 128, 1);
    broker = gpuBroker;
    nthreads += gpuBroker->GetNstream() + 1;
#else
    std::cerr << "Error: Coprocessor processing requested but support was not "
                 "enabled\n";
#endif
  }
  GeantPropagator *prop =
      GeantPropagator::Instance(ntotal, nbuffered, nthreads);
  if (broker)
    prop->SetTaskBroker(broker);
  // prop->fApplication = fApplication;
  prop->SetNminThreshold(5 * nthreads);
  prop->SetMonitored(GeantPropagator::kMonQueue, true & (!performance));
  prop->SetMonitored(GeantPropagator::kMonMemory, false & (!performance));
  prop->SetMonitored(GeantPropagator::kMonBasketsPerVol,
                     false & (!performance));
  prop->SetMonitored(GeantPropagator::kMonVectors, false & (!performance));
  prop->SetMonitored(GeantPropagator::kMonConcurrency, false & (!performance));
  prop->SetMonitored(GeantPropagator::kMonTracksPerEvent,
                     false & (!performance));
  prop->SetMonitored(GeantPropagator::kMonTracks, false & (!performance));
  bool graphics = (prop->GetMonFeatures()) ? true : false;
  prop->fUseMonitoring = graphics;
  // Value from individual vector
  prop->fPriorityThr = individual.GetPriority(individual);
  printf("Debugging RunCMS.C: priority value = %f\n", prop->fPriorityThr);
  // Value from individual vector
  prop->fNperBasket =
      individual.GetVector(individual); // Initial vector size (tunable)
  printf("Debugging RunCMS.C: vector value = %d\n", prop->fNperBasket);
  // Value from individual vector
  prop->fMaxPerBasket = 64; // Maximum vector size (tunable)
  prop->fMaxRes = 4000;
  if (performance)
    prop->fMaxRes = 0;
  prop->fEmin = 0.001; // [1 MeV] energy cut
  prop->fEmax = 0.01;  // 10 MeV
  prop->fProcess = new TTabPhysProcess("tab_phys", xsec, fstate);
  std::string s = "pp14TeVminbias.root";
  prop->fPrimaryGenerator = new HepMCGenerator(s);
  // Value from individual vector
  prop->fLearnSteps = individual.GetSteps(individual);
  printf("Debugging RunCMS.C: learning steps value = %d\n", prop->fLearnSteps);
  if (performance)
    prop->fLearnSteps = 0;
  CMSApplication *app = new CMSApplication();
  app->SetScoreType(CMSApplication::kScore);
  if (performance)
    app->SetScoreType(CMSApplication::kNoScore);
  prop->fApplication = app;
  prop->fDebugEvt = 0;
  prop->fDebugTrk = 0;
  prop->fDebugStp = 0;
  prop->fDebugRep = 10;
  prop->fUseStdScoring = true;
  if (performance)
    prop->fUseStdScoring = false;
  prop->fUseMonitoring = graphics;
  prop->PropagatorGeom(geomfile, nthreads, graphics);
#ifdef ENABLE_PERFMON
  perfcontrol.Stop();
#endif
  fitness->SetMemorySwitch(false);
  fitness->TemporarySolution();
  fitness->LogMemoryFitness("fitness.root");
  individual.SetFitness(0, prop->fTimer->RealTime());
  individual.SetFitness(1, -(prop->fNprimaries.load()));
  individual.SetFitness(2, fitness->GetmaxMemResident());
#ifdef ENABLE_PERFMON
  individual.SetFitness(3, perfcontrol.GetNInstructions());
  individual.SetFitness(4, perfcontrol.GetBranchMisses());
  perfcontrol.printSummary();
#endif
  delete prop;
  delete fitness;
  return individual.GetFitnessVector();
}

int main(int argc, char *argv[]) {
  // Function definition
  Functions *geantv = new Functions();
  geantv->SetIntervalGeantV();
  std::cout << "-==============================================-" << std::endl;
  geantv->PrintLimit(geantv->fInterval);
  std::cout << "-==============================================-" << std::endl;
  // Algorithm  definition
  AlgorithmNSGA *nsga2 = new AlgorithmNSGA();
  nsga2->SetPCross(0.5);
  nsga2->SetPMut(0.7);
  nsga2->SetGenTotalNumber(2);
  nsga2->SetNCons(0); // First version will be constrainless
  nsga2->SetNParam(6);
#ifdef ENABLE_PERFMON
  nsga2->SetNObjectives(5); // Memory, Time , etc.
#else
  nsga2->SetNObjectives(3);
#endif
  nsga2->SetCrowdingObj(false);
  nsga2->SetPopulationSize(4);
  nsga2->SetEtaMut(10);
  nsga2->SetEtaCross(10);
  nsga2->SetEpsilonC(0.7);
  nsga2->SetLimit(geantv->fInterval);
  nsga2->SetFunction(&CMSApp);
  nsga2->Initialize();
  nsga2->Evolution();
  return 0;
}
#else
int main(int argc, char *argv[]) {
  std::cout << "Geant-V based test: enable Geant-V in cmake flags and re-run "
               "compilation"
            << std::endl;
}
#endif
