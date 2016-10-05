#include "GATest.h"
#include "algorithms/GANSGA2.h"
#include "generic/Functions.h"
#include "generic/Population.h"
#include "generic/TGenes.h"
#include "instrumentation/GeantVFitness.h"
#include "output/HistogramManager.h"
#ifdef ENABLE_GEANTVVV
#include "problem/FastSimGeantV.h"

class Nsga2 : public GATest {
public:
};

TEST_F(Nsga2, SolvingGeantVFastSimProblem) {
  geantvmoop::FastSimGeantV runc;
  geantvmoop::GANSGA2<geantvmoop::FastSimGeantV> nsga(runc);
  nsga.fPopulationSize = 10;
  nsga.fMaxGeneration = 100;
  nsga.SolvePF();
}
#endif

