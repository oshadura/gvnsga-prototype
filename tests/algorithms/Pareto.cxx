#include "GATest.h"
#include "generic/PF.h"
#include "generic/Population.h"
#include "generic/TGenes.h"

#include "problem/DTLZ2.h"

class Pareto : public GATest {
public:
  geantvmoop::Population<geantvmoop::DTLZ2, 10> pop;
  geantvmoop::PF<geantvmoop::DTLZ2, 10> fFront;
};


TEST_F(Pareto, CheckPF) {
  auto fResult = geantvmoop::PF<geantvmoop::DTLZ2,10>::ParetoFrontND(pop);
  for (auto i : pop) {
    fFront.Add(i);
  }
  EXPECT_EQ(fResult, fFront.GetPopulation());
}

