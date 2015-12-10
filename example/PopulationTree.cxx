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

#include "Rtypes.h"
#include "TGeoManager.h"

#include "GunGenerator.h"
#include "HepMCGenerator.h"
#include "TaskBroker.h"
#include "WorkloadManager.h"
#include "GeantPropagator.h"
#include "TTabPhysProcess.h"
#include "CMSApplication.h"

#ifndef COPROCESSOR_REQUEST
#define COPROCESSOR_REQUEST false
#endif

void CMSApp(Genes<Double_t> *individual) {
  ///////////////////// Old part from counters experiments

  // const Events events {
  //  hwcounters::cache::L1::DCA,
  //  hwcounters::cache::L1::DCH,
  //  hwcounters::cache::L1::DCM };
  // auto counter = PerfStat(events);
  // counter.start();

  ///////////////////// Values from original function (runCMS.C)
  bool performance = true;
  const char *geomfile = "cms2015.root";
  const char *xsec = "xsec_FTFP_BERT_G496p02_1mev.root";
  const char *fstate = "fstate_FTFP_BERT_G496p02_1mev.root";
  bool coprocessor = COPROCESSOR_REQUEST;

  ///////////////////// Values should be taken from [Genes, map(Genes,Limits)]
  ///////////////////// Original macros
  //
  // Commenting line for compilation purposes, after we will get value from
  // individual vector
  // int nthreads = ncputhreads;
  int nthreads = 4;
  // Value from individual vector
  int ntotal = 10; // Number of events to be transported
  // Value from individual vector
  int nbuffered = 5; // Number of buffered events (tunable [1,ntotal])
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
  prop->fPriorityThr = 0.1;
  // Value from individual vector
  prop->fNperBasket = 16; // Initial vector size (tunable)
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
  prop->fLearnSteps = 100000;
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
  delete prop;
  // counter.stop;
  //////////////////// Write a Fitness [vector, map(FitnessValue,Constraint)]
}

int main() {
  // Function
  Functions *geantv = Functions::Instance();
  //geantv->SetInterval(); // don't work because we initialize fNparam after...
  geantv->fInterval.push_back(make_pair(1,10));
  geantv->fInterval.push_back(make_pair(1,10));
  std::vector<std::pair<Double_t,Double_t>> Interval = geantv->fInterval;
  std::cout << "Check what we create as a limit vector: [" << &Interval << "]"<< std::endl;
  // Algorithm  definition
  AlgorithmNSGA *nsga2 = new AlgorithmNSGA();
  nsga2->SetPCross(0.5);
  nsga2->SetEtaCross(0.7);
  nsga2->SetGenTotalNumber(5);
  nsga2->SetNCons(0); // First version will be constrainless
  nsga2->SetNParam(2); // blablabla - see Genes.h
  nsga2->SetNObjectives(2); // Memory, Time
  //nsga2->SetInterval(); // Testing intervals between [0,100]
  nsga2->SetPopulationSize(5);
  nsga2->SetPMut(0.7);
  nsga2->SetEtaMut(0.7);
  nsga2->SetEpsilonC(0.7);
  nsga2->SetFunction(&CMSApp);
  nsga2->Initialize();
  nsga2->Evolution();
  // Test population initialization
  /////////////////////////////////
  //Population<Double_t> fNsgaPop = new Population<Double_t>();
  //fNsgaPop.SetPopulationSize(20);
  //fNsgaPop.Build();
  //std::cout << fNsgaPop << std::endl;
  // fNsgaPop->WritePopulationTree(*fNsgaPop, "Population.root");
  // fNsgaPop->PrintTree("Population.root", "Population");
  /////////////////////////////////
  return 0;
}
