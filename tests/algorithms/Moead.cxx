#include "GATest.h"
#include "problem/DTLZ1.h"
#include "generic/Population.h"
#include "generic/Functions.h"
#include "generic/TGenes.h"
#include "output/HistogramManager.h"
#include "algorithms/GAMOEAD.h"
#include "instrumentation/GeantVFitness.h"
#include <boost/math/constants/constants.hpp>

class Moead : public GATest {
public:
};

// Wrong implementation of test..
TEST_F(Moead, SolvingProblem) {
	/**
  geantvmoop::DTLZ1 dtlz1;
  geantvmoop::GAMOEAD<geantvmoop::DTLZ1> moead(dtlz1);
  moead.fPopulationSize = 10;
  moead.SolvePF();
  */
}