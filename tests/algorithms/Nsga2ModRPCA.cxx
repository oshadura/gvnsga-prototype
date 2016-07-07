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
#include "algorithms/GANSGA2RPCA.h"
#include "instrumentation/GeantVFitness.h"

class Nsga2ModRPCA : public GATest {
public:
};

/*
TEST_F(Nsga2ModRPCA, SolvingDTLZ1Problem) {
  geantvmoop::DTLZ1 dtlz1;
  geantvmoop::GANSGA2ModRPCA<geantvmoop::DTLZ1> nsga2modrpca(dtlz1);
  nsga2modrpca.fPopulationSize = 100;
  nsga2modrpca.fMaxGeneration = 100;
  nsga2modrpca.SolvePF();
}
*/

/*
TEST_F(Nsga2ModRPCA, SolvingDTLZ2Problem) {
  geantvmoop::DTLZ2 dtlz2;
  geantvmoop::GANSGA2ModRPCA<geantvmoop::DTLZ2> nsga(dtlz2);
  nsga.fPopulationSize = 1000;
  nsga.fMaxGeneration = 100;
  nsga.SolvePF();
}
*/

/*
TEST_F(Nsga2ModRPCA, SolvingDTLZ3Problem) {
  geantvmoop::DTLZ3 dtlz3;
  geantvmoop::GANSGA2ModRPCA<geantvmoop::DTLZ3> nsga(dtlz3);
  nsga.fPopulationSize = 1000;
  nsga.fMaxGeneration = 100;
  nsga.SolvePF();
}

TEST_F(Nsga2ModRPCA, SolvingDTLZ4Problem) {
  geantvmoop::DTLZ4 dtlz4;
  geantvmoop::GANSGA2ModRPCA<geantvmoop::DTLZ4> nsga(dtlz4);
  nsga.fPopulationSize = 100;
  nsga.fMaxGeneration = 100;
  nsga.SolvePF();
}
*/
/*
TEST_F(Nsga2ModRPCA, SolvingDTLZ5Problem) {
  geantvmoop::DTLZ5 dtlz5;
  geantvmoop::GANSGA2ModRPCA<geantvmoop::DTLZ5> nsga(dtlz5);
  nsga.fPopulationSize = 10;
  nsga.fMaxGeneration = 10;
  nsga.SolvePF();
}

TEST_F(Nsga2ModRPCA, SolvingDTLZ6Problem) {
  geantvmoop::DTLZ6 dtlz6;
  geantvmoop::GANSGA2ModRPCA<geantvmoop::DTLZ6> nsga(dtlz6);
  nsga.fPopulationSize = 10;
  nsga.fMaxGeneration = 10;
  nsga.SolvePF();
}

TEST_F(Nsga2ModRPCA, SolvingDTLZ7Problem) {
  geantvmoop::DTLZ7 dtlz7;
  geantvmoop::GANSGA2ModRPCA<geantvmoop::DTLZ7> nsga(dtlz7);
  nsga.fPopulationSize = 10;
  nsga.fMaxGeneration = 10;
  nsga.SolvePF();
}
*/
