#include "GATest.h"
#include "PFMWatch.h"
#include "PAPIWatch.h"
#include "generic/Population.h"
#include "generic/Functions.h"
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
