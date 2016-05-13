#include "GATest.h"
#include "problem/ProblemDTLZ1.h"
#include "generic/Population.h"
#include "generic/Functions.h"
#include "generic/TGenes.h"
#include "output/HistogramManager.h"
#include "algorithms/MOEAD.h"
#include "instrumentation/GeantVFitness.h"
#include <boost/math/constants/constants.hpp>


class Moead : public GATest {
public:
	geantvmoop::ProblemDTLZ1 dtlz1;
};

/*

TEST_F(Moead, SolvingProblem) {
  geantvmoop::ProblemDTLZ1 dtlz1;
  geantvmoop::MOEAD<geantvmoop::ProblemDTLZ1> moead(dtlz1);
  moead.fPopulationSize = 1;
  moead.SolvePF();
  }
*/