#include "GATest.h"
#include "generic/Population.h"
#include "generic/TGenes.h"
#include "generic/PF.h"

#include "problem/DTLZ1.h"

class Pareto : public GATest {
public:
  geantvmoop::Population<geantvmoop::DTLZ1> pop{5};
  geantvmoop::PF<geantvmoop::DTLZ1> fFront;
};

TEST_F(Pareto, CheckPF) {
  auto fResult = geantvmoop::PF<geantvmoop::DTLZ1>::ParetoFrontND(pop);
  for (auto i : pop) {
    fFront.Add(i);
  }
  EXPECT_EQ(fResult, fFront.GetPopulation());
}