#include "GATest.h"
#include "algorithms/GANSGA2.h"
#include "algorithms/GANSGA2UPCA.h"
#include "generic/Functions.h"
#include "generic/Population.h"
#include "generic/TGenes.h"
#include "instrumentation/GeantVFitness.h"
#include "output/HistogramManager.h"
#include "problem/RunGeantV.h"

class Nsga2 : public GATest {
public:
};

#ifdef ENABLE_GEANTVVV
TEST_F(Nsga2, SolvingGeantVProblem) {
  geantvmoop::RunGeantV runc;
  geantvmoop::GANSGA2<geantvmoop::RunGeantV> nsga(runc);
  nsga.fPopulationSize = 100;
  nsga.fMaxGeneration = 100;
  nsga.SolvePF();
}
#endif

TEST_F(Nsga2, SolvingGeantVProblemUPCA) {
  geantvmoop::RunGeantV runc;
  geantvmoop::GANSGA2UPCA<geantvmoop::RunGeantV> nsgaupca(runc);
  nsgaupca.fPopulationSize = 50;
  nsgaupca.fMaxGeneration = 100;
  nsgaupca.SolvePF();
}
