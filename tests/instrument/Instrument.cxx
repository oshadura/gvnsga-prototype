#ifdef ENABLE_PAPI
#ifdef ENABLE_PERF

#include "GATest.h"
#include "PAPIWatch.h"
#include "PFMWatch.h"
#include "generic/Functions.h"
#include "generic/Population.h"
#include "generic/TGenes.h"
#include "problem/DTLZ2.h"

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

#endif
#endif
