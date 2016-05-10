#include "GATest.h"
#include "problem/ProblemDTLZ1.h"
#include "generic/Population.h"
#include "generic/Functions.h"
#include "output/HistogramManager.h"
#include "generic/TGenes.h"
#include "algorithms/NSGA.h"
#include "instrumentation/GeantVFitness.h"

#include <boost/math/constants/constants.hpp>


class Nsga2 : public GATest {
public:
	ProblemDTLZ1 dtlz1;
};

TEST_F(Nsga2, SolvingProblem) {
  NSGA<ProblemDTLZ1> nsga(dtlz1);
  nsga.fPopulationSize = 1;
  nsga.SolvePF();
  }


TEST_F(Nsga2, TestEvaluationOfOneIndividual){
}
