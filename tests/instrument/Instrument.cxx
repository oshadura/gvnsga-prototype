#include "GATest.h"
#include "PAPIWatch.h"
#include "PFMWatch.h"
#include "instrumentation/CPUManager.h"
//#include "Memory.h"
#include "generic/Functions.h"
#include "generic/Population.h"
#include "generic/TGenes.h"
#include "problem/DTLZ2.h"
// Externals pap-wrap from third-externals
#ifdef ENABLE_PAPI
#include "PapiCollectors.h"
#include "papi_wrap.h"
#endif

#include <iostream>
#include <chrono>
#include <thread>

/*
class Instrument : public GATest {
public:
#ifdef ENABLE_PAPI
  PAPIWatch papi;
#endif
  PFMWatch pfw;
  int i, j;
};

#ifdef ENABLE_PAPI

TEST_F(Instrument, CheckingPapi) {
  papi.setPapiEvents();
  papi.startPapi();
  static int x[4000][4000];
  for (i = 0; i < 4000; i++) {
    for (j = 0; j < 4000; j++) {
      x[j][i] = i + j;
    }
  }
  papi.stopPapi();
  papi.printPapiResults();
  papi.resetPapi();
}

#endif

TEST_F(Instrument, CheckingPerf) {
  pfw.Start();
  static int x[4000][4000];
  for (i = 0; i < 4000; i++) {
    for (j = 0; j < 4000; j++) {
      x[j][i] = i + j;
    }
  }
  pfw.Stop();
  pfw.printSummary();
}

#ifdef ENABLE_PAPI

TEST_F(Instrument, CheckingPapiWrap) {
  int handle = pw_new_collector("PapiTestOnPopulation");
  pw_start_collector(handle);
  static int x[4000][4000];
  for (i = 0; i < 4000; i++) {
    for (j = 0; j < 4000; j++) {
      x[j][i] = i + j;
    }
  }
  pw_stop_collector(handle);
  pw_print();
  pw_print_table();
}

#endif

TEST_F(Instrument, CheckMemory) {
  size_t currentSize = getCurrentRSS();
  size_t peakSize = getPeakRSS();
  std::cout << "Current RSS Size: " << currentSize << std::endl;
  std::cout << "Peak RSS Size: " << peakSize << std::endl;
}

TEST_F(Instrument, CheckCPU) {
  CPUManager cpumgr;
  cpumgr.InitCPU();
  hwloc_topology_t topology;
  double nbcores, ccores;
  hwloc_topology_init(&topology); // initialization
  hwloc_topology_load(topology);  // actual detection
  nbcores = hwloc_get_nbobjs_by_type(topology, HWLOC_OBJ_PU);
  hwloc_topology_destroy(topology);
  static int x[4000][4000];
  for (i = 0; i < 4000; i++) {
    for (j = 0; j < 4000; j++) {
      x[j][i] = i + j;
    }
  }
  ccores = nbcores - cpumgr.GetCurrentValueCPU() / 100 * nbcores; // just a test
std::cout << " Number of total free cores " << ccores << std::endl;
}

*/