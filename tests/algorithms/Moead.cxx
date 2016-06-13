#include "GATest.h"
#include "problem/DTLZ1.h"
#include "problem/DTLZ2.h"
#include "problem/DTLZ3.h"
#include "problem/DTLZ4.h"
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
*/

/*
TEST_F(Moead, SolvingProblemDTLZ2) {
  geantvmoop::DTLZ2 dtlz2;
  geantvmoop::GAMOEAD<geantvmoop::DTLZ2> moead(dtlz2);
  moead.fPopulationSize = 100;
  moead.fMaxGeneration = 100;
  moead.SolvePF();
  moead.PrintImpl(std::cout);
}

/*
TEST_F(Moead, SolvingProblemDTLZ4) {
  geantvmoop::DTLZ4 dtlz4;
  geantvmoop::GAMOEAD<geantvmoop::DTLZ4> moead(dtlz4);
  moead.fPopulationSize = 100;
  moead.fMaxGeneration = 100;
  moead.SolvePF();
  moead.PrintImpl(std::cout);
}
*/

