#include "GATest.h"
#include "problem/DTLZ1.h"
#include "problem/DTLZ2.h"
#include "problem/DTLZ3.h"
#include "problem/DTLZ4.h"
#include "problem/DTLZ5.h"
#include "problem/DTLZ6.h"
#include "problem/DTLZ7.h"
#include "generic/Population.h"
#include "generic/Functions.h"
#include "generic/TGenes.h"
#include "output/HistogramManager.h"
#include "algorithms/GANSGA2UPCA.h"
#include "instrumentation/GeantVFitness.h"

class Nsga2upca : public GATest {
public:
};

/*
TEST_F(Nsga2upca, SolvingDTLZ1Problem) {
  geantvmoop::DTLZ1 dtlz1;
  geantvmoop::GANSGA2UPCA<geantvmoop::DTLZ1> nsga2upca(dtlz1);
  nsga2upca.fPopulationSize = 1000;
  nsga2upca.fMaxGeneration = 40;
  nsga2upca.SolvePF();
}

/*
TEST_F(Nsga2upca, SolvingDTLZ2Problem) {
  geantvmoop::DTLZ2 dtlz2;
  geantvmoop::GANSGA2UPCA<geantvmoop::DTLZ2> nsga(dtlz2);
  nsga.fPopulationSize = 1000;
  nsga.fMaxGeneration = 100;
  nsga.SolvePF();
}
*/

/*
TEST_F(Nsga2upca, SolvingDTLZ3Problem) {
  geantvmoop::DTLZ3 dtlz3;
  geantvmoop::GANSGA2UPCA<geantvmoop::DTLZ3> nsga(dtlz3);
  nsga.fPopulationSize = 1000;
  nsga.fMaxGeneration = 100;
  nsga.SolvePF();
}
*/


/*
TEST_F(Nsga2upca, SolvingDTLZ4Problem) {
  geantvmoop::DTLZ4 dtlz4;
  geantvmoop::GANSGA2upca<geantvmoop::DTLZ4> nsga(dtlz4);
  nsga.fPopulationSize = 100;
  nsga.fMaxGeneration = 100;
  nsga.SolvePF();
}
*/
/*
TEST_F(Nsga2upca, SolvingDTLZ5Problem) {
  geantvmoop::DTLZ5 dtlz5;
  geantvmoop::GANSGA2upca<geantvmoop::DTLZ5> nsga(dtlz5);
  nsga.fPopulationSize = 10;
  nsga.fMaxGeneration = 10;
  nsga.SolvePF();
}

TEST_F(Nsga2upca, SolvingDTLZ6Problem) {
  geantvmoop::DTLZ6 dtlz6;
  geantvmoop::GANSGA2upca<geantvmoop::DTLZ6> nsga(dtlz6);
  nsga.fPopulationSize = 10;
  nsga.fMaxGeneration = 10;
  nsga.SolvePF();
}

TEST_F(Nsga2upca, SolvingDTLZ7Problem) {
  geantvmoop::DTLZ7 dtlz7;
  geantvmoop::GANSGA2upca<geantvmoop::DTLZ7> nsga(dtlz7);
  nsga.fPopulationSize = 10;
  nsga.fMaxGeneration = 10;
  nsga.SolvePF();
}
*/
