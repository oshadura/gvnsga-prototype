#include "GATest.h"
#include "problem/DTLZ1.h"
#include "generic/Population.h"
#include "generic/Functions.h"
#include "generic/TGenes.h"
#include "output/HistogramManager.h"
#include "algorithms/GANSGA2.h"
#include "instrumentation/GeantVFitness.h"

class Nsga2 : public GATest {
public:
};

TEST_F(Nsga2, SolvingDTLZ1Problem) {
  geantvmoop::DTLZ1 dtlz1;
  geantvmoop::GANSGA2<geantvmoop::DTLZ1> nsga(dtlz1);
  nsga.fPopulationSize = 10;
  nsga.fMaxGeneration = 10;
  nsga.SolvePF();
}

