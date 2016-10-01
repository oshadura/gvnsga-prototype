#include "GATest.h"
#include "algorithms/GANSGA2.h"
#include "generic/Functions.h"
#include "generic/Population.h"
#include "generic/TGenes.h"
#include "instrumentation/GeantVFitness.h"
#include "output/HistogramManager.h"
#include "problem/CMSGeantV.h"

class Nsga2 : public GATest {
public:
};

#ifdef ENABLE_GEANTV
TEST_F(Nsga2, SolvingGeantVCMSProblem) {
  geantvmoop::CMSGeantV runc;
  geantvmoop::GANSGA2<geantvmoop::CMSGeantV> nsga(runc);
  nsga.fPopulationSize = 10;
  nsga.fMaxGeneration = 100;
  nsga.SolvePF();
}
#endif

