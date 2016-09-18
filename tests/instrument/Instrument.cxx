#include "GATest.h"
#include "PAPIWatch.h"
#include "PFMWatch.h"
#include "Memory.h"
#include "generic/Functions.h"
#include "generic/Population.h"
#include "generic/TGenes.h"
#include "problem/DTLZ2.h"
// Externals pap-wrap from third-externals
#include "PapiCollectors.h"
#include "papi_wrap.h"

class Instrument : public GATest {
public:
  PAPIWatch papi;
  PFMWatch pfw;
};

TEST_F(Instrument, CheckingPapi) {
  papi.setPapiEvents();
  papi.startPapi();
  geantvmoop::TGenes<geantvmoop::DTLZ2> i;
  papi.stopPapi();
  papi.printPapiResults();
}

TEST_F(Instrument, CheckingPerf) {
  pfw.Start();
  geantvmoop::TGenes<geantvmoop::DTLZ2> i;
  pfw.Stop();
  pfw.printSummary();
}

TEST_F(Instrument, CheckingPapiWrap){
  int handle = pw_new_collector("PapiTestOnPopulation");
  pw_start_collector(handle);
  geantvmoop::TGenes<geantvmoop::DTLZ2> i;
  pw_stop_collector(handle);
  pw_print();
  pw_print_table();
}

TEST_F(Instrument, CheckMemory){
  size_t currentSize = getCurrentRSS( );
  size_t peakSize    = getPeakRSS( );
  std::cout << "Current RSS Size: " << currentSize << std::endl;
  std::cout << "Peak RSS Size: " << peakSize << std::endl;
}