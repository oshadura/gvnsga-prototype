#include "GATest.h"
#include "problem/DTLZ2.h"
#include "gaoperators/GACrossover.h"
#include "gaoperators/GASBXCrossover.h"

class CrossoverSBXTest : public GATest {
public:
  typedef geantvmoop::GAVector<geantvmoop::GADouble, 12> SBXInput;

  geantvmoop::individual_t<geantvmoop::DTLZ2> fIndividual1 =
      std::make_shared<geantvmoop::TGenes<geantvmoop::DTLZ2>>(SBXInput{12, geantvmoop::GADouble(0, 0, 10000)});
  geantvmoop::individual_t<geantvmoop::DTLZ2> fIndividual2 =
      std::make_shared<geantvmoop::TGenes<geantvmoop::DTLZ2>>(SBXInput{12, geantvmoop::GADouble(0, 0, 10000)});
};
