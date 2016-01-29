#include <cmath>
#include <iostream>  // std::cout
#include <iterator>  // std::ostream_iterator
#include <vector>    // std::vector
#include <algorithm> // std::copy

#include "Population.h"
#include "Functions.h"
#include "HistogramManager.h"
#include "TGenes.h"
#include "AlgorithmNSGA.h"
#include "GeantVFitness.h"

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

const char *geomfile = "ExN03.root";
const char *xsec = "xsec_FTFP_BERT.root";
const char *fstate = "fstate_FTFP_BERT.root";
bool performance = true;
bool coprocessor = COPROCESSOR_REQUEST;
int nthreads = 4;
int ntotal = 10;
int nbuffered = 5;
GeantPropagator *prop =
      GeantPropagator::Instance(ntotal, nbuffered, nthreads);



void runAppLoop(Genes<Double_t> &individual) {
  prop->fProcess = new TTabPhysProcess("tab_phys", xsec, fstate);
  // Should be generated each time because of prop.Clear()
  prop->fApplication = new ExN03Application();
  //////////////////////////////////////////////////
  prop->fNthreads = individual.GetThread(individual);
  printf("Debugging Run.C: thread value = %d\n", prop->fNthreads);
  Int_t nthread = prop->fNtotal =
      individual.GetAllev(individual); // Number of events to be transported
  printf("Debugging Run.C: all events value = %d\n", prop->fNtotal);
  prop->fNevents = individual.GetBuffev(
      individual); // Number of buffered events (tunable [1,ntotal])
  printf("Debugging Run.C: buffered particles value = %d\n", prop->fNevents);
  prop->fNaverage = 500; 
  prop->fPriorityThr = individual.GetPriority(individual);
  printf("Debugging Run.C: priority value = %f\n", prop->fPriorityThr);
  prop->fNperBasket =
      individual.GetVector(individual); // Initial vector size (tunable)
  printf("Debugging Run.C: vector value = %d\n", prop->fNperBasket);
  prop->fMaxPerBasket = individual.GetMaxVector(individual); // Maximum vector size (tunable)
  printf("Debugging Run.C: max vector value = %d\n", prop->fMaxPerBasket);
  // for vector physics -OFF now
  // prop->fVectorPhysicsProcess = new GVectorPhysicsProcess(prop->fEmin,
  // nthreads);
  prop->fLearnSteps = individual.GetSteps(individual);
  printf("Debugging Run.C: learning steps value = %d\n", prop->fLearnSteps);
  if (performance)
    prop->fLearnSteps = 0;
  prop->fPrimaryGenerator =
      new GunGenerator(prop->fNaverage, 11, prop->fEmax, -8, 0, 0, 1, 0, 0);
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
  prop->fUseMonitoring = false;
  prop->fUseAppMonitoring = false;
  prop->PropagatorGeom(geomfile, nthread, false);
  individual.SetFitness(0, prop->fTimer->RealTime());
  // fitness.HistOutputFitness();
  prop->Clear();
  return;
}

int main(int argc, char *argv[]) {
printf("First initialized value in Run.C: thread value = %d\n", nthreads);
printf("First initialized value in Run.C: ntotal = %d\n", ntotal);
printf("First initialized value in Run.C: nbuffered = %d\n", nbuffered);
#ifdef ENABLE_PERFMON
  PFMWatch perfcontrol;
#endif
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
  if (broker)
    prop->SetTaskBroker(broker);
  // Monitor different features
  prop->SetNminThreshold(5 * nthreads);
  prop->SetMonitored(GeantPropagator::kMonQueue, true & (!performance));
  prop->SetMonitored(GeantPropagator::kMonMemory, false & (!performance));
  prop->SetMonitored(GeantPropagator::kMonBasketsPerVol,
                     false & (!performance));
  prop->SetMonitored(GeantPropagator::kMonVectors, false & (!performance));
  prop->SetMonitored(GeantPropagator::kMonConcurrency, false & (!performance));
  prop->SetMonitored(GeantPropagator::kMonTracksPerEvent,
                     false & (!performance));
  prop->fNaverage = 500; // Average number of tracks per event
  // Threshold for prioritizing events (tunable [0, 1], normally <0.1)
  // If set to 0 takes the default value of 0.01
  prop->fPriorityThr = 0.01;
  printf("First initialized value in Run.C: priority value = %f\n", prop->fPriorityThr);
  // Initial vector size, this is no longer an important model parameter,
  // because is gets dynamically modified to accomodate the track flow
  prop->fNperBasket = 64; // Initial vector size (tunable)
  printf("First initialized value in Run.C: fNperBasket = %d\n", prop->fNperBasket);
  // This is now the most important parameter for memory considerations
  prop->fMaxPerBasket = 256; // Maximum vector size (tunable)
  prop->fEmin = 3.E-6;       // [3 KeV] energy cut
  prop->fEmax = 0.03; // [30MeV] used for now to select particle gun energy
  // for vector physics -OFF now
  // prop->fVectorPhysicsProcess = new GVectorPhysicsProcess(prop->fEmin,
  // nthreads);
  prop->fPrimaryGenerator =
      new GunGenerator(prop->fNaverage, 11, prop->fEmax, -8, 0, 0, 1, 0, 0);
  // Number of steps for learning phase (tunable [0, 1e6])
  // if set to 0 disable learning phase
  prop->fLearnSteps = 1000;
  printf("irst initialized value in Run.C: learning steps value = %d\n", prop->fLearnSteps);
  if (performance)
    prop->fLearnSteps = 0;
  /////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////
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
  nsga2->SetGenTotalNumber(40);
  nsga2->SetNCons(0); // First version will be constrainless
  nsga2->SetNParam(7);
  nsga2->SetNObjectives(2); // Memory, Time
  // nsga2->SetInterval(); // Testing intervals between [0,100]
  nsga2->SetCrowdingObj(false);
  nsga2->SetPopulationSize(4);
  nsga2->SetEtaMut(10);
  nsga2->SetEtaCross(10);
  nsga2->SetEpsilonC(0.7);
  nsga2->SetLimit(geantv->fInterval);
  // nsga2->SetFunction(&CMSApp);
  nsga2->SetFunction(&runAppLoop);
  nsga2->Initialize();
  nsga2->Evolution();
  return 0;
}
