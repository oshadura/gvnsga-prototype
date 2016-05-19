#include "GATest.h"
#include "problem/RunGeantV.h"
#include "generic/Population.h"
#include "generic/Functions.h"
#include "generic/TGenes.h"
#include "output/HistogramManager.h"
#include "algorithms/GANSGA2.h"
#include "instrumentation/GeantVFitness.h"
#include <boost/math/constants/constants.hpp>

class Nsga2 : public GATest {
public:
};

TEST_F(Nsga2, SolvingGeantVProblem) {
  geantvmoop::RunGeantV runc;
  geantvmoop::GANSGA2<geantvmoop::RunGeantV> nsga(runc);
  nsga.fPopulationSize = 10;
  nsga.SolvePF();
}
