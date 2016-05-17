#include "GATest.h"
#include "problem/DTLZ1.h"
#include "gaoperators/GACrossover.h"
#include "gaoperators/GASBXCrossover.h"

class CrossoverSBXTest : public GATest {
public:
  typedef geantvmoop::GAVector<geantvmoop::GADouble> SBXInput;
  geantvmoop::individual_t<geantvmoop::DTLZ1> fIndividual1 =
      std::make_shared<geantvmoop::TGenes<geantvmoop::DTLZ1>>(
          SBXInput{3, geantvmoop::GADouble(0, 0, 10000)});
  geantvmoop::individual_t<geantvmoop::DTLZ1> fIndividual2 =
      std::make_shared<geantvmoop::TGenes<geantvmoop::DTLZ1>>(
          SBXInput{3, geantvmoop::GADouble(0, 0, 10000)});
};

/*
TEST_F(CrossoverSBXTest, SBX){
        ASSERT_NO_THROW(geantvmoop::SBXCrossover::CrossoverGA(fIndividual1,fIndividual2));
}
*/
