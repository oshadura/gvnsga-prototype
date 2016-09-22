#include "gtest/gtest.h"
#include "GATest.h"
#include "generic/Population.h"
#include "problem/DTLZ2.h"
#include "output/HistogramManager.h"

class HistoOutput : public GATest {
public:
  geantvmoop::DTLZ2 dtlz2;
  geantvmoop::Population<geantvmoop::DTLZ2> pop{5};
};

/*
TEST_F(HistoOutput, SimpleHistoOutputofPopulation) {
	std::cout << pop;
  bool IsTrue = HistogramManager<geantvmoop::DTLZ2>::GetInstance().HistoFill(
      pop, "testpopulation.root", 0);
  ASSERT_EQ(IsTrue, true);
}
*/
