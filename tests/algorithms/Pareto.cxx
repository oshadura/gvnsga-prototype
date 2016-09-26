#include "GATest.h"
#include "generic/Population.h"
#include "generic/TGenes.h"
#include "generic/PF.h"

#include "problem/DTLZ2.h"

class Pareto : public GATest {
public:
  geantvmoop::Population<geantvmoop::DTLZ2> pop{5};
  geantvmoop::PF<geantvmoop::DTLZ2> fFront;
};

/*
TEST_F(Pareto, CheckPF) {
  auto fResult = geantvmoop::PF<geantvmoop::DTLZ2>::ParetoFrontND(pop);
  for (auto i : pop) {
    fFront.Add(i);
  }
  EXPECT_EQ(fResult, fFront.GetPopulation());
}
*/