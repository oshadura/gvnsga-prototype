#include "Population.h"
#include "Functions.h"
#include "HistogramManager.h"
#include "Genes.h"
#include "AlgorithmNSGA.h"


#include "GeantPropagator.h"
#include "GeantScheduler.h"
#include "WorkloadManager.h"
#include "PhysicsProcess.h"
#include "GeantVApplication.h"
#include "ExN03Application.h"
#include "GeantFactoryStore.h"
#include "GeantTrack.h"
#include "GunGenerator.h"
#include "TTabPhysProcess.h"
#include "HepMCGenerator.h"


#ifndef COPROCESSOR_REQUEST
#define COPROCESSOR_REQUEST false
#endif

void Run(int nthreads=4,
         bool performance=true,
        const char *geomfile="ExN03.root",
        const char *xsec="xsec_FTFP_BERT.root",
        const char *fstate="fstate_FTFP_BERT.root",
         bool coprocessor = COPROCESSOR_REQUEST){  
  //const Events events {
  //  hwcounters::cache::L1::DCA,
  //  hwcounters::cache::L1::DCH,
  //  hwcounters::cache::L1::DCM };
  //auto counter = PerfStat(events);
  //counter.start();
  int ntotal   = 50;  // Number of events to be transported
  int nbuffered  = 10;   // Number of buffered events (tunable [1,ntotal])
  TGeoManager::Import(geomfile); 
  GeantPropagator *prop = GeantPropagator::Instance(ntotal, nbuffered, nthreads);
   // Monitor different features
  prop->SetNminThreshold(5*nthreads);
  if (coprocessor) {
#ifdef GEANTCUDA_REPLACE
      CoprocessorBroker *gpuBroker = new CoprocessorBroker();
      gpuBroker->CudaSetup(32,128,1);
      prop->SetTaskBroker(gpuBroker);
#else
      std::cerr << "Error: Coprocessor processing requested but support was not enabled\n";
#endif
   }
  prop->SetMonitored(GeantPropagator::kMonQueue,          true & (!performance));
  prop->SetMonitored(GeantPropagator::kMonMemory,         false & (!performance));
  prop->SetMonitored(GeantPropagator::kMonBasketsPerVol,  false & (!performance));
  prop->SetMonitored(GeantPropagator::kMonVectors,        false & (!performance));
  prop->SetMonitored(GeantPropagator::kMonConcurrency,    false & (!performance));
  prop->SetMonitored(GeantPropagator::kMonTracksPerEvent, false & (!performance));
  bool graphics = (prop->GetMonFeatures()) ? true : false;
  prop->fUseMonitoring = graphics;
  prop->fNaverage = 500;   // Average number of tracks per event
   // Threshold for prioritizing events (tunable [0, 1], normally <0.1)
   // If set to 0 takes the default value of 0.01
  prop->fPriorityThr = 0.05;
   // Initial vector size, this is no longer an important model parameter, 
   // because is gets dynamically modified to accomodate the track flow
  prop->fNperBasket = 16;   // Initial vector size (tunable)
   // This is now the most important parameter for memory considerations
  prop->fMaxPerBasket = 256;   // Maximum vector size (tunable)
  prop->fEmin = 3.E-6; // [3 KeV] energy cut
  prop->fEmax = 0.03;  // [30MeV] used for now to select particle gun energy
   // Create the tab. phys process.
  prop->fProcess = new TTabPhysProcess("tab_phys", xsec, fstate);
   // for vector physics -OFF now
   // prop->fVectorPhysicsProcess = new GVectorPhysicsProcess(prop->fEmin, nthreads);
  prop->fPrimaryGenerator = new GunGenerator(prop->fNaverage, 11, prop->fEmax, -8, 0, 0, 1, 0, 0);
   // Number of steps for learning phase (tunable [0, 1e6])
   // if set to 0 disable learning phase
  prop->fLearnSteps = 0;
  if (performance) prop->fLearnSteps = 0;
  prop->fApplication = new ExN03Application();
   // Activate I/O
  prop->fFillTree = false;
// Activate debugging using -DBUG_HUNT=ON in your cmake build
  prop->fDebugEvt = 0;
  prop->fDebugTrk = 0;
  prop->fDebugStp = 0;
  prop->fDebugRep = 10;
// Activate standard scoring   
  prop->fUseStdScoring = true;
  if (performance) prop->fUseStdScoring = false;
   // Monitor the application
  prop->fUseAppMonitoring = false;
  prop->PropagatorGeom(geomfile, nthreads, graphics);
  delete prop;
  //counter.stop;
}

int main() {
  // Function
  Functions *geantv = Functions::Instance();
  geantv->SetNCons(0);
  geantv->SetNParams(6);

  // Setup limits for generation parameters

  // Setup constrain = 0
  //
  //Setup Function
  //geantv->SetFunction(*(Run));
  //geantv->SetFunctionOpt(*(RunC),4);

  // Setup genes generation

  // Algorithm  definition
  AlgorithmNSGA *nsga2 = AlgorithmNSGA::Instance();
  nsga2->SetPopulationSize(5);
  nsga2->SetPCross(0.5);
  nsga2->SetEtaCross(0.7);
  nsga2->SetEtaMut(0.3);
  nsga2->SetPMut(0.5);
  printf("Population size (presetuped) = %d\n",nsga2->GetPopulationSetupSize());
  printf("Probability of crossover = %g\n",nsga2->GetPCross());
  printf ("Eta values for crossover (crossover rate) = %g\n",nsga2->GetEtaCross());
  printf ("Probability  for mutation = %g\n",nsga2->GetPMut());
  printf ("Eta values for mutation (mutation rate)= %g\n",nsga2->GetEtaMut());
  // Setup population initialization
  Population *fNsgaPop;
  //Testing function
  //fNsgaPop->WritePopulationTree(*fNsgaPop, "Population");
  //fNsgaPop->PrintTree("Population", "Population");
  fNsgaPop->Build();
  // Algorithm run
  nsga2->Initialize();
  nsga2->Evolution();

  // Missing logging process...
  return 0;
}
