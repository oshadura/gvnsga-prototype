#include "GATest.h"
#include "algorithms/GANSGA2UPCA.h"
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
#include "problem/Fonseca.h"
#include "problem/Kursawe.h"
#include "problem/Shwefel.h"
#include "problem/UFC9.h"

class Nsga2upca : public GATest {
public:
};
/*
TEST_F(Nsga2upca, SolvingKursaweProblem) {
  geantvmoop::Kursawe f;
  geantvmoop::GANSGA2UPCA<geantvmoop::Kursawe> nsga2upca(f);
  nsga2upca.fPopulationSize = 1000;
  nsga2upca.fMaxGeneration = 500;
  nsga2upca.SolvePF();
}
TEST_F(Nsga2upca, SolvingDTLZ1Problem) {
  geantvmoop::DTLZ1 dtlz1;
  geantvmoop::GANSGA2UPCA<geantvmoop::DTLZ1> nsga(dtlz1);
  nsga.fPopulationSize = 50;
  nsga.fMaxGeneration = 100;
  nsga.SolvePF();
}
TEST_F(Nsga2upca, SolvingDTLZ3Problem) {
  geantvmoop::DTLZ3 dtlz3;
  geantvmoop::GANSGA2UPCA<geantvmoop::DTLZ3> nsga(dtlz3);
  nsga.fPopulationSize = 1000;
  nsga.fMaxGeneration = 100;
  nsga.SolvePF();
}

TEST_F(Nsga2upca, SolvingShwefelProblem) {
  geantvmoop::Shwefel shw;
  geantvmoop::GANSGA2UPCA<geantvmoop::Shwefel> nsga(shw);
  nsga.fPopulationSize = 1000;
  nsga.fMaxGeneration = 100;
  nsga.SolvePF();
}
TEST_F(Nsga2upca, SolvingDTLZ4Problem) {
  geantvmoop::DTLZ4 dtlz4;
  geantvmoop::GANSGA2UPCA<geantvmoop::DTLZ4> nsga(dtlz4);
  nsga.fPopulationSize = 1000;
  nsga.fMaxGeneration = 100;
  nsga.SolvePF();
}

TEST_F(Nsga2upca, SolvingDTLZ5Problem) {
  geantvmoop::DTLZ5 dtlz5;
  geantvmoop::GANSGA2UPCA<geantvmoop::DTLZ5> nsga(dtlz5);
  nsga.fPopulationSize = 100;
  nsga.fMaxGeneration = 200;
  nsga.SolvePF();
}

TEST_F(Nsga2upca, SolvingDTLZ6Problem) {
  geantvmoop::DTLZ6 dtlz6;
  geantvmoop::GANSGA2UPCA<geantvmoop::DTLZ6> nsga(dtlz6);
  nsga.fPopulationSize = 1000;
  nsga.fMaxGeneration = 100;
  nsga.SolvePF();
}

TEST_F(Nsga2upca, SolvingDTLZ7Problem) {
  geantvmoop::DTLZ7 dtlz7;
  geantvmoop::GANSGA2UPCA<geantvmoop::DTLZ7> nsga(dtlz7);
  nsga.fPopulationSize = 1000;
  nsga.fMaxGeneration = 100;
  nsga.SolvePF();
}

/*
TEST_F(Nsga2upca, SolvingDTLZ7Problem) {
  geantvmoop::UFC9 dtlz7;
  geantvmoop::GANSGA2UPCA<geantvmoop::UFC9> nsga(dtlz7);
  nsga.fPopulationSize = 1000;
  nsga.fMaxGeneration = 100;
  nsga.SolvePF();
}
*/
