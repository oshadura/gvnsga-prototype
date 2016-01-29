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
#include "PFMWatch.h"

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

bool performance = true;
const char *geomfile = "cms2015.root";
const char *xsec = "xsec_FTFP_BERT_G496p02_1mev.root";
const char *fstate = "fstate_FTFP_BERT_G496p02_1mev.root";
bool coprocessor = COPROCESSOR_REQUEST;
int nthreads = 4;
int ntotal = 10;
int nbuffered = 5;
GeantPropagator *prop = GeantPropagator::Instance(ntotal, nbuffered, nthreads);


void CMSApp(Genes<Double_t> &individual) {
  //GeantVFitness fitness;
  std::cout << "Lets pass it to GeantV propagator.." << std::endl;
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
  // Value from individual vector
  prop->fPriorityThr = individual.GetPriority(individual);
  printf("Debugging RunCMS.C: priority value = %f\n", prop->fPriorityThr);
  // Value from individual vector
  prop->fNperBasket =
      individual.GetVector(individual); // Initial vector size (tunable)
  printf("Debugging RunCMS.C: vector value = %d\n", prop->fNperBasket);
  // Value from individual vector
  prop->fMaxPerBasket = individual.GetMaxVector(individual);; // Maximum vector size (tunable)
  printf("Debugging RunCMS.C: vector value = %d\n", prop->fMaxPerBasket);
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
  prop->PropagatorGeom(geomfile, nthreads, prop->fUseMonitoring);
  delete prop;
  individual.SetFitness(0, prop->fTimer->RealTime());
  //fitness.HistOutputFitness();
  return;
}

int main(int argc, char *argv[]) {
  printf("First initialized value in Run.C: thread value = %d\n", nthreads);
  printf("First initialized value in Run.C: ntotal = %d\n", ntotal);
  printf("First initialized value in Run.C: nbuffered = %d\n", nbuffered);

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
  prop->fUseMonitoring = (prop->GetMonFeatures()) ? true : false;
  prop->fMaxRes = 4000;
  if (performance)
    prop->fMaxRes = 0;
  prop->fEmin = 0.001; // [1 MeV] energy cut
  prop->fEmax = 0.01;  // 10 MeV
  prop->fProcess = new TTabPhysProcess("tab_phys", xsec, fstate);
  std::string s = "pp14TeVminbias.root";
  prop->fPrimaryGenerator = new HepMCGenerator(s);
  ///////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
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
  nsga2->SetGenTotalNumber(3);
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
  nsga2->SetFunction(&CMSApp);
  nsga2->Initialize();
  nsga2->Evolution();
  return 0;
}
