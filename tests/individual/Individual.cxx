#include "GATest.h"
#include "problem/DTLZ1.h"
#include "generic/Population.h"
#include "generic/Functions.h"
#include "generic/TGenes.h"

class Individual : public GATest {
public:
};

TEST_F(Individual, GeneratingIndividual) {
	geantvmoop::TGenes<geantvmoop::DTLZ1> i;
}
