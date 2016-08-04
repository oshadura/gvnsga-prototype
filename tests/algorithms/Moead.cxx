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
#include "algorithms/GAMOEAD.h"
#include "instrumentation/GeantVFitness.h"
#include <boost/math/constants/constants.hpp>

class Moead : public GATest {
public:
};

/*
TEST_F(Moead, SolvingProblemDTLZ1) {
  geantvmoop::DTLZ1 dtlz1;
  geantvmoop::GAMOEAD<geantvmoop::DTLZ1> moead(dtlz1);
  moead.fPopulationSize = 100;
  moead.fMaxGeneration = 100;
  moead.SolvePF();
  moead.PrintImpl(std::cout);
}

TEST_F(Moead, SolvingProblemDTLZ2) {
  geantvmoop::DTLZ2 dtlz2;
  geantvmoop::GAMOEAD<geantvmoop::DTLZ2> moead(dtlz2);
  moead.fPopulationSize = 100;
  moead.fMaxGeneration = 100;
  moead.SolvePF();
  moead.PrintImpl(std::cout);
}

TEST_F(Moead, SolvingProblemDTLZ3) {
  geantvmoop::DTLZ3 dtlz3;
  geantvmoop::GAMOEAD<geantvmoop::DTLZ3> moead(dtlz3);
  moead.fPopulationSize = 100;
  moead.fMaxGeneration = 100;
  moead.SolvePF();
  moead.PrintImpl(std::cout);
}

TEST_F(Moead, SolvingProblemDTLZ4) {
  geantvmoop::DTLZ4 dtlz4;
  geantvmoop::GAMOEAD<geantvmoop::DTLZ4> moead(dtlz4);
  moead.fPopulationSize = 100;
  moead.fMaxGeneration = 100;
  moead.SolvePF();
  moead.PrintImpl(std::cout);
}

TEST_F(Moead, SolvingProblemDTLZ5) {
  geantvmoop::DTLZ5 dtlz5;
  geantvmoop::GAMOEAD<geantvmoop::DTLZ5> moead(dtlz5);
  moead.fPopulationSize = 100;
  moead.fMaxGeneration = 100;
  moead.SolvePF();
  moead.PrintImpl(std::cout);
}

TEST_F(Moead, SolvingProblemDTLZ6) {
  geantvmoop::DTLZ6 dtlz6;
  geantvmoop::GAMOEAD<geantvmoop::DTLZ6> moead(dtlz6);
  moead.fPopulationSize = 100;
  moead.fMaxGeneration = 100;
  moead.SolvePF();
  moead.PrintImpl(std::cout);
}

TEST_F(Moead, SolvingProblemDTLZ7) {
  geantvmoop::DTLZ7 dtlz7;
  geantvmoop::GAMOEAD<geantvmoop::DTLZ7> moead(dtlz7);
  moead.fPopulationSize = 100;
  moead.fMaxGeneration = 100;
  moead.SolvePF();
  moead.PrintImpl(std::cout);
}
*/