#include "GATest.h"
#include "problem/DTLZ2.h"
#include "generic/Population.h"
#include "generic/Functions.h"
#include "generic/TGenes.h"

class Individual : public GATest {
public:
};

TEST_F(Individual, GeneratingIndividual) {
  geantvmoop::TGenes<geantvmoop::DTLZ2> i;
}


/*
TEST_F(Individual, GeneratingIndividualDifferent) {
	        typename F::Input tmpinput;
        typename F::Output tmpoutput;
  auto forkedindividual = std::make_shared<TGenes<F> >(tmpinput, tmpoutput);
  auto indv = (*forkedindividual).GetInput();
  std::cout << "--------------------------------------" << std::endl;
}
*/

TEST_F(Individual, GenerationPopulation) {
  geantvmoop::Population<geantvmoop::DTLZ2> pop{ 5 };
}
