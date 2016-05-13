#include "GATest.h"
#include "problem/ProblemDTLZ1.h"
#include "gaoperators/Crossover.h"
#include "gaoperators/SBXCrossover.h"


class CrossoverSBXTest : public GATest {
public:
	typedef geantvmoop::GaVector<geantvmoop::RandomDouble> SBXInput;
	geantvmoop::individual_t<geantvmoop::ProblemDTLZ1> fIndividual1 = std::make_shared<geantvmoop::Genes<geantvmoop::ProblemDTLZ1>>(SBXInput{3, geantvmoop::RandomDouble(0,0,10000)});
	geantvmoop::individual_t<geantvmoop::ProblemDTLZ1> fIndividual2 = std::make_shared<geantvmoop::Genes<geantvmoop::ProblemDTLZ1>>(SBXInput{3, geantvmoop::RandomDouble(0,0,10000)});
};

/*
TEST_F(CrossoverSBXTest, SBX){
	ASSERT_NO_THROW(geantvmoop::SBXCrossover::CrossoverGA(fIndividual1,fIndividual2));
}
*/
