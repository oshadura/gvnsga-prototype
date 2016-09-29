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
#include "algorithms/GANSGA2.h"
#include "instrumentation/GeantVFitness.h"

class Nsga2 : public GATest {
public:
};

TEST_F(Nsga2, SolvingDTLZ1Problem) {
  geantvmoop::DTLZ1 dtlz1;
  geantvmoop::GANSGA2<geantvmoop::DTLZ1, 100, 7> nsga(dtlz1);
  nsga.fMaxGeneration = 500;
  nsga.SolvePF();
}

TEST_F(Nsga2, SolvingDTLZ2Problem) {
  geantvmoop::DTLZ2 dtlz2;
  geantvmoop::GANSGA2<geantvmoop::DTLZ2, 100, 12> nsga(dtlz2);
  nsga.fMaxGeneration = 500;
  nsga.SolvePF();
}

TEST_F(Nsga2, SolvingDTLZ3Problem) {
  geantvmoop::DTLZ3 dtlz3;
  geantvmoop::GANSGA2<geantvmoop::DTLZ3, 100, 7> nsga(dtlz3);
  nsga.fMaxGeneration = 500;
  nsga.SolvePF();
}

TEST_F(Nsga2, SolvingDTLZ4Problem) {
  geantvmoop::DTLZ4 dtlz4;
  geantvmoop::GANSGA2<geantvmoop::DTLZ4, 100, 12> nsga(dtlz4);
  nsga.fMaxGeneration = 500;
  nsga.SolvePF();
}

TEST_F(Nsga2, SolvingDTLZ5Problem) {
  geantvmoop::DTLZ5 dtlz5;
  geantvmoop::GANSGA2<geantvmoop::DTLZ5, 100, 7> nsga(dtlz5);
  nsga.fMaxGeneration = 500;
  nsga.SolvePF();
}

TEST_F(Nsga2, SolvingDTLZ6Problem) {
  geantvmoop::DTLZ6 dtlz6;
  geantvmoop::GANSGA2<geantvmoop::DTLZ6, 100, 7> nsga(dtlz6);
  nsga.fMaxGeneration = 500;
  nsga.SolvePF();
}

TEST_F(Nsga2, SolvingDTLZ7Problem) {
  geantvmoop::DTLZ7 dtlz7;
  geantvmoop::GANSGA2<geantvmoop::DTLZ7, 100, 12> nsga(dtlz7);
  nsga.fMaxGeneration = 500;
  nsga.SolvePF();
}


