#include "GATest.h"
#include "algorithms/GANSGA3.h"
#include "generic/Functions.h"
#include "generic/Population.h"
#include "generic/TGenes.h"
#include "instrumentation/GeantVFitness.h"
#include "output/HistogramManager.h"
#include "problem/DTLZ1.h"
#include "problem/DTLZ2.h"
#include "problem/DTLZ3.h"
#include "problem/DTLZ4.h"
#include "problem/DTLZ5.h"
#include "problem/DTLZ6.h"
#include "problem/DTLZ7.h"

class Nsga3 : public GATest {
public:
};
/*
TEST_F(NSGA3, SolvingDTLZ1Problem) {
  geantvmoop::DTLZ1 dtlz1;
  geantvmoop::GANSGA3<geantvmoop::DTLZ1> nsga(dtlz1);
  nsga.fPopulationSize = 1000;
  nsga.fMaxGeneration = 100;
  nsga.SolvePF();
}

TEST_F(NSGA3, SolvingDTLZ2Problem) {
  geantvmoop::DTLZ2 dtlz2;
  geantvmoop::GANSGA3<geantvmoop::DTLZ2> nsga(dtlz2);
  nsga.fPopulationSize = 1000;
  nsga.fMaxGeneration = 100;
  nsga.SolvePF();
}

TEST_F(NSGA3, SolvingDTLZ3Problem) {
  geantvmoop::DTLZ3 dtlz3;
  geantvmoop::GANSGA3<geantvmoop::DTLZ3> nsga(dtlz3);
  nsga.fPopulationSize = 1000;
  nsga.fMaxGeneration = 100;
  nsga.SolvePF();
}

TEST_F(NSGA3, SolvingDTLZ4Problem) {
  geantvmoop::DTLZ4 dtlz4;
  geantvmoop::GANSGA3<geantvmoop::DTLZ4> nsga(dtlz4);
  nsga.fPopulationSize = 1000;
  nsga.fMaxGeneration = 100;
  nsga.SolvePF();
}

TEST_F(NSGA3, SolvingDTLZ5Problem) {
  geantvmoop::DTLZ5 dtlz5;
  geantvmoop::GANSGA3<geantvmoop::DTLZ5> nsga(dtlz5);
  nsga.fPopulationSize = 1000;
  nsga.fMaxGeneration = 100;
  nsga.SolvePF();
}

TEST_F(NSGA3, SolvingDTLZ6Problem) {
  geantvmoop::DTLZ6 dtlz6;
  geantvmoop::GANSGA3<geantvmoop::DTLZ6> nsga(dtlz6);
  nsga.fPopulationSize = 1000;
  nsga.fMaxGeneration = 100;
  nsga.SolvePF();
}

TEST_F(NSGA3, SolvingDTLZ7Problem) {
  geantvmoop::DTLZ7 dtlz7;
  geantvmoop::GANSGA3<geantvmoop::DTLZ7> nsga(dtlz7);
  nsga.fPopulationSize = 1000;
  nsga.fMaxGeneration = 100;
  nsga.SolvePF();
}
*/
