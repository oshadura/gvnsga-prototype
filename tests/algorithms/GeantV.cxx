#include "GATest.h"
#include "algorithms/GANSGA2.h"
#include "generic/Functions.h"
#include "generic/Population.h"
#include "generic/TGenes.h"
#include "instrumentation/GeantVFitness.h"
#include "output/HistogramManager.h"
#include "problem/RunGeantV.h"

class Nsga2 : public GATest {
public:
};

#ifdef ENABLE_GEANTV


TEST_F(Nsga2, SolvingGeantVProblem) {
  geantvmoop::RunGeantV runc;
  geantvmoop::GANSGA2<geantvmoop::RunGeantV> nsga(runc);
  nsga.fPopulationSize = 1;
  nsga.SolvePF();
}


#endif
