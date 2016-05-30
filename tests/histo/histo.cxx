#include "gtest/gtest.h"
#include "GATest.h"
#include "generic/Population.h"
#include "problem/DTLZ1.h"
#include "output/HistogramManager.h"

class HistoOutput : public GATest {
public:
  geantvmoop::DTLZ1 dtlz1;
  geantvmoop::Population<geantvmoop::DTLZ1> pop{5};
};


TEST_F(HistoOutput, SimpleHistoOutputofPopulation) {
  geantvmoop::Population<geantvmoop::DTLZ1> population;
  bool IsTrue = HistogramManager<geantvmoop::DTLZ1>::GetInstance().HistoFill(population, "testpopulation.root");
  ASSERT_EQ(IsTrue, true);
}